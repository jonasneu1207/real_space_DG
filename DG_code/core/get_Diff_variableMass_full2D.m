function [A, rhs, info] = get_Diff_variableMass_full2D(mat, p, boundary)
%GET_DIFF_VARIABLEMASS_FULL2D Full-2D kinetic operator for m=m(X,Y).
%
% The constant-mass Wigner equation contains only the mixed derivatives
% d_X*d_rho_x and d_Y*d_rho_y. That simplification is not valid across a
% material interface. For the Hermitian BenDaniel-Duke Hamiltonian
%
%   T = -(hbar^2/2) div( m(r)^(-1) grad )
%
% the density-matrix kinetic term for one Cartesian direction is evaluated
% before any Wigner/Fourier transformation as
%
%   K_m f = -1i*hbar/2 * [ D_+ (a_+ D_+ f)
%                         - D_- (a_- D_- f) ],
%
%   D_+ = 0.5*d_R + d_rho,   D_- = 0.5*d_R - d_rho,
%   a_+ = 1/m(R+rho/2),      a_- = 1/m(R-rho/2).
%
% This file applies that expression independently in x and y. d_R is the
% same central rectangular-DG derivative used by get_Diff_full2D, and
% d_rho is the existing cell-centred FV derivative. The implementation is
% therefore a direct discretization of the variable-mass von-Neumann
% operator, rather than a local rescaling of the constant-mass transport.
%
% To retain the established Source/Drain and Y-wall numerical fluxes, the
% operator is written as
%
%   K_m = K_ref + (K_m - K_ref).
%
% K_ref is built by get_Diff_full2D with scalar reference masses. The
% correction below contains the complete endpoint-mass expression and is
% exactly zero for a spatially constant mass equal to the reference mass.
% Consequently the old solver is recovered to roundoff in that limit.
% Open-boundary data still use the homogeneous-contact reference operator;
% this is the standard lead assumption and requires the mass to have
% reached its contact value at the outer X faces.
%
% Global DOF order is unchanged:
%   f(iRhoX,iRhoY,iX,iY), with rho_x fastest.

if nargin < 3 || isempty(boundary)
    boundary = get_Boundary_full2D(mat, p, p.EfL, p.EfR, []);
end

params = getDGParams(mat);
referenceMassX = referenceMassRelative(mat, params, 'X');
referenceMassY = referenceMassRelative(mat, params, 'Y');

% Build the already validated flux/boundary operator with explicitly known
% reference masses. Only the new correction depends on the full mass maps.
referenceMat = mat;
referenceMat.me_x_ch = referenceMassX;
referenceMat.me_y_ch = referenceMassY;
[AReference, rhs, referenceInfo] = get_Diff_full2D( ...
    referenceMat, p, boundary);

data = variableMassData(mat, p, params, referenceInfo, ...
    referenceMassX, referenceMassY);

if referenceInfo.fullMatrixAssembled
    A = sparse(AReference + assembleCorrection(data));
else
    A = [];
    rhs = [];
end

% Preserve the public interface of get_Diff_full2D so the system builder,
% Krylov solvers and all existing preconditioners can use either path.
info = referenceInfo;
info.massModel = 'position-dependent-bdd';
info.massDiscretization = 'endpoint BenDaniel-Duke von-Neumann operator';
info.referenceTransport = referenceInfo;
info.referenceMassRelative = [referenceMassX, referenceMassY];
info.mass = data.massInfo;
info.reason = variableMassReason(referenceInfo, p.index.nTotal);
info.apply = @(u) applyVariableMass(referenceInfo, data, u);
info.makeGPUApply = @() makeGPUVariableMassApply(referenceInfo, data);
info.assemble = @() sparse(referenceInfo.assemble() + assembleCorrection(data));
info.assembleRhs = @() referenceInfo.assembleRhs();
info.getDiagonal = @() referenceInfo.getDiagonal() ...
    + getCorrectionDiagonal(data);
info.getRowAbsSum = @() referenceInfo.getRowAbsSum() ...
    + getCorrectionRowAbsBound(data);
info.getRelativeBlockData = @() makeVariableMassBlockData( ...
    referenceInfo, data);
end

function data = variableMassData(mat, p, params, referenceInfo, massX, massY)
constants = physicalConstants;

data = struct;
data.nCenter = p.dg.nCenterDof;
data.nRelative = p.relative.nDof;
data.nTotal = p.index.nTotal;
data.centerDerivativeX = centerDerivative(p, 'X');
data.centerDerivativeY = centerDerivative(p, 'Y');
data.relativeDerivativeX = kron(speye(p.relative.NrhoY), ...
    p.relative.DrhoX);
data.relativeDerivativeY = kron(p.relative.DrhoY, ...
    speye(p.relative.NrhoX));

% The sampled coefficients are dimensionless inverse relative masses
% b=m0/m. This avoids storing numbers of order 1e31 and improves numerical
% cancellation in homogeneous regions. Since K_ref already contains b_ref,
% only b_+ - b_ref and b_- - b_ref are retained.
[data.massXPlus, data.massXMinus, massInfoX] = ...
    endpointInverseMassCorrection(mat, p, params, 'X', massX);
[data.massYPlus, data.massYMinus, massInfoY] = ...
    endpointInverseMassCorrection(mat, p, params, 'Y', massY);

% If Qx=hbar/(m0*m_ref), this coefficient is -i*hbar/(2*m0).
% Expressing it through the reference scale also keeps the existing
% full2D_diff_scale and optional Q overrides effective on the new path.
data.coefficientX = -0.5i*referenceInfo.transportScales.Qx*massX;
data.coefficientY = -0.5i*referenceInfo.transportScales.Qy*massY;
data.hasXCorrection = nnz(data.massXPlus) > 0 ...
    || nnz(data.massXMinus) > 0;
data.hasYCorrection = nnz(data.massYPlus) > 0 ...
    || nnz(data.massYMinus) > 0;

data.massInfo = struct;
data.massInfo.X = massInfoX;
data.massInfo.Y = massInfoY;
data.massInfo.referenceRelative = [massX, massY];
data.massInfo.endpointDefinition = ...
    'b_+/- = m0/m(R +/- rho/2), sampled from mat.me_x and mat.me_y';
data.massInfo.operator = ...
    '-i*hbar/2 [D_+(m_+^{-1}D_+) - D_-(m_-^{-1}D_-)]';
data.massInfo.constantMassRegression = [ ...
    'The correction is exactly zero when mat.me_x/me_y equal the ', ...
    'selected scalar reference masses.'];
data.massInfo.contactAssumption = [ ...
    'Source/Drain characteristic data use homogeneous lead reference ', ...
    'masses; mass profiles should be constant in the contact extensions.'];
data.massInfo.storageBytes = 8*(numel(data.massXPlus) ...
    + numel(data.massXMinus) + numel(data.massYPlus) ...
    + numel(data.massYMinus));
data.constants = constants;
end

function C = centerDerivative(p, axisName)
%CENTERDERIVATIVE Reproduce the central first-derivative DG operator.
% This is intentionally algebraically identical to the H1 construction in
% get_Diff_full2D so the constant-mass regression is exact.
switch upper(axisName)
    case 'X'
        H = oneDimensionalCentralDerivative(p.dg.X);
        C = kron(speye(p.dg.Y.nDof), H);
    case 'Y'
        H = oneDimensionalCentralDerivative(p.dg.Y);
        C = kron(H, speye(p.dg.X.nDof));
end
C = sparse(C);
end

function H = oneDimensionalCentralDerivative(axis)
nLocal = axis.nLocal;
nElement = axis.nElements;
invM = sparse(inv(full(axis.M)));
P1 = sparse(1, 1, 1, nLocal, nLocal);
P3 = sparse(nLocal, nLocal, 1, nLocal, nLocal);
P4 = sparse(nLocal, 1, 1, nLocal, nLocal);
P2 = sparse(1, nLocal, 1, nLocal, nLocal);

K1 = invM*(axis.S + 0.5*P1 - 0.5*P3);
K3 = -0.5*invM*P2;
K5 = 0.5*invM*P4;
upperScale = spdiags(1./axis.jacobian(:), 1, nElement, nElement);
lowerScale = spdiags(1./axis.jacobian(2:end).', -1, ...
    nElement, nElement);
diagScale = spdiags(1./axis.jacobian(:), 0, nElement, nElement);
H = kron(upperScale, K5) + kron(lowerScale, K3) ...
    + kron(diagScale, K1);
H = sparse(H);
end

function [plusCorrection, minusCorrection, info] = ...
        endpointInverseMassCorrection(mat, p, params, axisName, referenceMass)
[inverseMass, gridX, gridY, fieldInfo] = inverseRelativeMassField( ...
    mat, p, axisName, referenceMass);
massInterpolant = griddedInterpolant({gridX, gridY}, inverseMass, ...
    'linear', 'nearest');

centerX = repmat(p.dg.X.nodes(:), p.dg.Y.nDof, 1);
centerY = kron(p.dg.Y.nodes(:), ones(p.dg.X.nDof, 1));
relativeX = repmat(p.relative.rhoX.cells(:), p.relative.NrhoY, 1);
relativeY = kron(p.relative.rhoY.cells(:), ...
    ones(p.relative.NrhoX, 1));

nRelative = p.relative.nDof;
nCenter = p.dg.nCenterDof;
plusCorrection = zeros(nRelative, nCenter);
minusCorrection = zeros(nRelative, nCenter);
targetEntries = max(1, floor(readParam(params, ...
    'full2D_variableMassChunkTargetEntries', 2e6)));
centerChunk = max(1, floor(targetEntries/nRelative));

for first = 1:centerChunk:nCenter
    last = min(nCenter, first + centerChunk - 1);
    ids = first:last;
    X = centerX(ids).';
    Y = centerY(ids).';
    plusCorrection(:, ids) = massInterpolant( ...
        X + 0.5*relativeX, Y + 0.5*relativeY) - 1/referenceMass;
    minusCorrection(:, ids) = massInterpolant( ...
        X - 0.5*relativeX, Y - 0.5*relativeY) - 1/referenceMass;
end

coefficientTolerance = readParam(params, ...
    'full2D_variableMassCoefficientTol', 1e-12);
absoluteTolerance = coefficientTolerance*max(1, 1/referenceMass);
plusCorrection(abs(plusCorrection) <= absoluteTolerance) = 0;
minusCorrection(abs(minusCorrection) <= absoluteTolerance) = 0;

leftIds = sub2ind([p.dg.X.nDof, p.dg.Y.nDof], ...
    ones(p.dg.Y.nDof, 1), (1:p.dg.Y.nDof).');
rightIds = sub2ind([p.dg.X.nDof, p.dg.Y.nDof], ...
    p.dg.X.nDof*ones(p.dg.Y.nDof, 1), (1:p.dg.Y.nDof).');

info = fieldInfo;
info.referenceRelativeMass = referenceMass;
info.inverseReferenceMass = 1/referenceMass;
info.plusCorrectionRange = rangeOf(plusCorrection);
info.minusCorrectionRange = rangeOf(minusCorrection);
info.maxOuterXFaceCorrection = max(abs([ ...
    plusCorrection(:, leftIds), minusCorrection(:, leftIds), ...
    plusCorrection(:, rightIds), minusCorrection(:, rightIds)]), [], 'all');
info.outerXFacesMatchReference = info.maxOuterXFaceCorrection ...
    <= absoluteTolerance;
info.interpolation = 'linear inverse mass, nearest contact extension';
info.chunkCenterDof = centerChunk;
end

function [inverseMass, gridX, gridY, info] = ...
        inverseRelativeMassField(mat, p, axisName, defaultMass)
fieldName = ['me_', lower(axisName)];
massField = readField(mat, fieldName, defaultMass);
massField = massFieldToXY(massField, mat, defaultMass);

constants = physicalConstants;
kgMask = abs(massField) < 1e-25;
massField(kgMask) = massField(kgMask)/constants.m0;
if any(~isfinite(massField), 'all') || any(massField <= 0, 'all')
    error('DG:Full2D:InvalidVariableMass', ...
        '%s must contain finite, strictly positive effective masses.', fieldName);
end
inverseMass = 1./massField;

gridX = coordinateVector(mat, 'x', 'dx', 'Nx')*p.coordinateScale;
gridY = coordinateVector(mat, 'y', 'dy', 'Ny')*p.coordinateScale;
[gridX, ix] = sort(gridX(:));
[gridY, iy] = sort(gridY(:));
inverseMass = full(inverseMass(ix, iy));

info = struct;
info.field = fieldName;
info.relativeMassRange = [min(massField, [], 'all'), ...
    max(massField, [], 'all')];
info.inverseRelativeMassRange = [min(inverseMass, [], 'all'), ...
    max(inverseMass, [], 'all')];
info.gridSize = size(inverseMass);
end

function field = massFieldToXY(field, mat, defaultValue)
if isempty(field)
    field = defaultValue;
end
if isscalar(field)
    field = repmat(field, mat.Nx, mat.Ny);
    return
end

if ~ismatrix(field)
    fieldSize = size(field);
    if fieldSize(end-1) ~= mat.Nx || fieldSize(end) ~= mat.Ny
        error('DG:Full2D:InvalidVariableMassSize', ...
            'Mass array must end in dimensions Nx-by-Ny.');
    end
    field = reshape(field, [], mat.Nx, mat.Ny);
    % The present Full-2D solver advances one valley. This matches the
    % existing contact-mode and transport paths, which use valley 1.
    field = squeeze(field(1, :, :));
else
    field = squeeze(field);
end

if isvector(field)
    if numel(field) == mat.Nx
        field = repmat(field(:), 1, mat.Ny);
    elseif numel(field) == mat.Ny
        field = repmat(field(:).', mat.Nx, 1);
    elseif numel(field) == mat.Nx*mat.Ny
        field = reshape(field, mat.Nx, mat.Ny);
    else
        error('DG:Full2D:InvalidVariableMassSize', ...
            'Mass vector has %d entries.', numel(field));
    end
elseif isequal(size(field), [mat.Nx, mat.Ny])
    % Already in the material-grid convention field(X,Y).
elseif isequal(size(field), [mat.Ny, mat.Nx])
    field = field.';
else
    error('DG:Full2D:InvalidVariableMassSize', ...
        'Mass field has size [%d,%d], expected [%d,%d].', ...
        size(field, 1), size(field, 2), mat.Nx, mat.Ny);
end
end

function y = applyVariableMass(referenceInfo, data, u)
if numel(u) ~= data.nTotal
    error('DG:Full2D:InvalidOperatorInput', ...
        'Input has %d entries, expected %d.', numel(u), data.nTotal);
end
Y = reshape(referenceInfo.apply(u), data.nRelative, data.nCenter);
U = reshape(u, data.nRelative, data.nCenter);
if data.hasXCorrection
    Y = Y + applyAxisCorrection(U, data.centerDerivativeX, ...
        data.relativeDerivativeX, data.massXPlus, data.massXMinus, ...
        data.coefficientX);
end
if data.hasYCorrection
    Y = Y + applyAxisCorrection(U, data.centerDerivativeY, ...
        data.relativeDerivativeY, data.massYPlus, data.massYMinus, ...
        data.coefficientY);
end
y = Y(:);
end

function Y = applyAxisCorrection(U, C, R, massPlus, massMinus, coefficient)
plusGradient = 0.5*applyCenter(C, U) + R*U;
minusGradient = 0.5*applyCenter(C, U) - R*U;
plusFlux = massPlus.*plusGradient;
minusFlux = massMinus.*minusGradient;
plusDivergence = 0.5*applyCenter(C, plusFlux) + R*plusFlux;
minusDivergence = 0.5*applyCenter(C, minusFlux) - R*minusFlux;
Y = coefficient*(plusDivergence - minusDivergence);
end

function Y = applyCenter(C, U)
if isa(U, 'gpuArray')
    Y = (C*U.').';
else
    Y = U*C.';
end
end

function applyGPU = makeGPUVariableMassApply(referenceInfo, data)
referenceApply = referenceInfo.makeGPUApply();
gpuData = data;
gpuData.centerDerivativeX = gpuArray(data.centerDerivativeX);
gpuData.centerDerivativeY = gpuArray(data.centerDerivativeY);
gpuData.relativeDerivativeX = gpuArray(data.relativeDerivativeX);
gpuData.relativeDerivativeY = gpuArray(data.relativeDerivativeY);
gpuData.massXPlus = gpuArray(data.massXPlus);
gpuData.massXMinus = gpuArray(data.massXMinus);
gpuData.massYPlus = gpuArray(data.massYPlus);
gpuData.massYMinus = gpuArray(data.massYMinus);
gpuReferenceInfo = struct('apply', referenceApply);
applyGPU = @(u) applyVariableMass(gpuReferenceInfo, gpuData, u);
end

function correction = assembleCorrection(data)
correction = spalloc(data.nTotal, data.nTotal, 0);
if data.hasXCorrection
    correction = correction + assembleAxisCorrection( ...
        data.centerDerivativeX, data.relativeDerivativeX, ...
        data.massXPlus, data.massXMinus, data.coefficientX, ...
        data.nCenter, data.nRelative);
end
if data.hasYCorrection
    correction = correction + assembleAxisCorrection( ...
        data.centerDerivativeY, data.relativeDerivativeY, ...
        data.massYPlus, data.massYMinus, data.coefficientY, ...
        data.nCenter, data.nRelative);
end
correction = sparse(correction);
end

function A = assembleAxisCorrection(C, R, massPlus, massMinus, ...
        coefficient, nCenter, nRelative)
Icenter = speye(nCenter);
Irelative = speye(nRelative);
centerPart = 0.5*kron(C, Irelative);
relativePart = kron(Icenter, R);
Dplus = centerPart + relativePart;
Dminus = centerPart - relativePart;
A = coefficient*(Dplus*spdiags(massPlus(:), 0, ...
    nCenter*nRelative, nCenter*nRelative)*Dplus ...
    - Dminus*spdiags(massMinus(:), 0, ...
    nCenter*nRelative, nCenter*nRelative)*Dminus);
A = sparse(A);
end

function diagonal = getCorrectionDiagonal(data)
diagonalMatrix = complex(zeros(data.nRelative, data.nCenter));
if data.hasXCorrection
    diagonalMatrix = diagonalMatrix + axisCorrectionDiagonal( ...
        data.centerDerivativeX, data.relativeDerivativeX, ...
        data.massXPlus, data.massXMinus, data.coefficientX);
end
if data.hasYCorrection
    diagonalMatrix = diagonalMatrix + axisCorrectionDiagonal( ...
        data.centerDerivativeY, data.relativeDerivativeY, ...
        data.massYPlus, data.massYMinus, data.coefficientY);
end
diagonal = diagonalMatrix(:);
end

function diagonal = axisCorrectionDiagonal(C, R, plusMass, minusMass, coefficient)
centerKernel = C.*C.';
relativeKernel = R.*R.';
cross = diag(R)*diag(C).';
plusDiagonal = 0.25*(plusMass*centerKernel.') ...
    + cross.*plusMass + relativeKernel*plusMass;
minusDiagonal = 0.25*(minusMass*centerKernel.') ...
    - cross.*minusMass + relativeKernel*minusMass;
diagonal = coefficient*(plusDiagonal - minusDiagonal);
end

function rowAbs = getCorrectionRowAbsBound(data)
rowAbsMatrix = zeros(data.nRelative, data.nCenter);
if data.hasXCorrection
    rowAbsMatrix = rowAbsMatrix + axisCorrectionRowAbsBound( ...
        data.centerDerivativeX, data.relativeDerivativeX, ...
        data.massXPlus, data.massXMinus, data.coefficientX);
end
if data.hasYCorrection
    rowAbsMatrix = rowAbsMatrix + axisCorrectionRowAbsBound( ...
        data.centerDerivativeY, data.relativeDerivativeY, ...
        data.massYPlus, data.massYMinus, data.coefficientY);
end
rowAbs = rowAbsMatrix(:);
end

function bound = axisCorrectionRowAbsBound(C, R, plusMass, minusMass, coefficient)
absC = abs(C);
absR = abs(R);
rowMagnitude = full(sum(absR, 2)) ...
    + 0.5*full(sum(absC, 2)).';
plusWeighted = abs(plusMass).*rowMagnitude;
minusWeighted = abs(minusMass).*rowMagnitude;
plusBound = absR*plusWeighted + 0.5*plusWeighted*absC.';
minusBound = absR*minusWeighted + 0.5*minusWeighted*absC.';
bound = abs(coefficient)*(plusBound + minusBound);
end

function blockData = makeVariableMassBlockData(referenceInfo, data)
referenceBlocks = referenceInfo.getRelativeBlockData();
blockData = struct;
blockData.nCenter = data.nCenter;
blockData.nRelative = data.nRelative;
blockData.nTotal = data.nTotal;
blockData.blockType = 'variable-mass-relative-rho-blocks-per-center-dof';
blockData.dofOrder = {'rho_x', 'rho_y', 'X-DG-DOF', 'Y-DG-DOF'};
blockData.note = ['Each block contains the established reference transport ', ...
    'plus the exact center-block diagonal of the endpoint-mass correction.'];
blockData.getBlock = @(centerId) variableMassBlock( ...
    referenceBlocks, data, centerId);
end

function block = variableMassBlock(referenceBlocks, data, centerId)
block = referenceBlocks.getBlock(centerId);
if data.hasXCorrection
    block = block + axisCorrectionBlock(data.centerDerivativeX, ...
        data.relativeDerivativeX, data.massXPlus, data.massXMinus, ...
        data.coefficientX, centerId);
end
if data.hasYCorrection
    block = block + axisCorrectionBlock(data.centerDerivativeY, ...
        data.relativeDerivativeY, data.massYPlus, data.massYMinus, ...
        data.coefficientY, centerId);
end
block = sparse(block);
end

function block = axisCorrectionBlock(C, R, plusMass, minusMass, ...
        coefficient, centerId)
plusBlock = derivativeCenterBlock(C, R, plusMass, centerId, 1);
minusBlock = derivativeCenterBlock(C, R, minusMass, centerId, -1);
block = coefficient*(plusBlock - minusBlock);
end

function block = derivativeCenterBlock(C, R, mass, centerId, signR)
massHere = mass(:, centerId);
M = spdiags(massHere, 0, size(R, 1), size(R, 2));
centerRoundTrip = full(C(centerId, :).*C(:, centerId).');
centerMass = mass*centerRoundTrip.';
block = R*M*R ...
    + 0.5*signR*C(centerId, centerId)*(M*R + R*M) ...
    + 0.25*spdiags(centerMass, 0, size(R, 1), size(R, 2));
block = sparse(block);
end

function mass = referenceMassRelative(mat, params, axisName)
fieldName = ['full2D_referenceMass', upper(axisName)];
defaultName = ['me_', lower(axisName), '_ch'];
defaultMass = readField(mat, defaultName, 0.041);
mass = readParam(params, fieldName, defaultMass);
if ~isscalar(mass) || ~isfinite(mass) || mass <= 0
    error('DG:Full2D:InvalidReferenceMass', ...
        '%s must be a finite positive scalar.', fieldName);
end
constants = physicalConstants;
if abs(mass) < 1e-25
    mass = mass/constants.m0;
end
end

function text = variableMassReason(referenceInfo, nTotal)
if referenceInfo.fullMatrixAssembled
    text = sprintf(['Variable-mass Full-2D kinetic matrix assembled for ', ...
        '%d DOFs using endpoint BenDaniel-Duke coefficients.'], nTotal);
else
    text = sprintf(['Variable-mass Full-2D kinetic operator kept ', ...
        'matrix-free for %d DOFs. Reference flux path: %s'], ...
        nTotal, referenceInfo.reason);
end
end

function limits = rangeOf(values)
limits = [min(values, [], 'all'), max(values, [], 'all')];
end

function coord = coordinateVector(mat, coordName, spacingName, countName)
if isfield(mat, coordName) && ~isempty(mat.(coordName))
    coord = mat.(coordName)(:);
else
    coord = (0:mat.(countName)-1).'*mat.(spacingName);
end
end

function params = getDGParams(mat)
params = struct;
if isfield(mat, 'dg') && isfield(mat.dg, 'params')
    params = mat.dg.params;
end
end

function value = readParam(params, name, defaultValue)
value = defaultValue;
if isstruct(params) && isfield(params, name) && ~isempty(params.(name))
    value = params.(name);
end
end

function value = readField(s, name, defaultValue)
value = defaultValue;
if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
    value = s.(name);
end
end

function c = physicalConstants
c.hbar = 1.054571817e-34;
c.m0 = 9.1093837015e-31;
end
