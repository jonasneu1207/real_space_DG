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
info.massDiscretization = ...
    'quadrature-projected endpoint BenDaniel-Duke von-Neumann operator';
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
[centerMass, centerMassInverse] = exactCenterMass(p);
data.centerDerivativeX = centerDerivative(p, 'X');
data.centerDerivativeY = centerDerivative(p, 'Y');

% The outer BDD divergence must be the adjoint in the same exact DG inner
% product that defines p.dg.X.M and p.dg.Y.M. In particular, the dense
% element mass matrix must not be replaced by diagonal LGL weights here.
% Doing that would mix a mass-lumped variable coefficient with the exact
% mass matrix used by the reference transport operator.
data.centerAdjointX = sparse(centerMassInverse ...
    * data.centerDerivativeX' * centerMass);
data.centerAdjointY = sparse(centerMassInverse ...
    * data.centerDerivativeY' * centerMass);
data.relativeDerivativeX = kron(speye(p.relative.NrhoY), ...
    p.relative.DrhoX);
data.relativeDerivativeY = kron(p.relative.DrhoY, ...
    speye(p.relative.NrhoX));
data.relativeAdjointX = data.relativeDerivativeX';
data.relativeAdjointY = data.relativeDerivativeY';

projection = massProjectionData(p, params);
projectionNames = fieldnames(projection);
for id = 1:numel(projectionNames)
    data.(projectionNames{id}) = projection.(projectionNames{id});
end

% The sampled coefficients are dimensionless inverse relative masses
% b=m0/m. This avoids storing numbers of order 1e31 and improves numerical
% cancellation in homogeneous regions. Since K_ref already contains b_ref,
% only b_+ - b_ref and b_- - b_ref are retained.
[data.massXPlus, data.massXMinus, massInfoX] = ...
    endpointInverseMassCorrection(mat, p, params, projection, 'X', massX);
[data.massYPlus, data.massYMinus, massInfoY] = ...
    endpointInverseMassCorrection(mat, p, params, projection, 'Y', massY);

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
    ['b_+/- = m0/m(R +/- rho/2), evaluated at tensor-product ', ...
     'center-coordinate Gauss-Lobatto points'];
data.massInfo.operator = ...
    '-i*hbar/2 [-D_+^*(m_+^{-1}D_+) + D_-^*(m_-^{-1}D_-)]';
data.massInfo.discreteDivergence = [ ...
    'The inverse mass is projected locally as M^(-1) E'' W diag(b) E. ', ...
    'The outer derivative is the exact DG/FV weighted adjoint of the ', ...
    'inner gradient, using the same non-diagonal element mass matrix as ', ...
    'the reference DG transport operator.'];
data.massInfo.projection = struct( ...
    'type', 'element-local strong DG multiplication', ...
    'formula', 'M^(-1) E'' W diag(b) E', ...
    'quadratureOrder', projection.quadratureOrder, ...
    'nCenterQuadrature', projection.nCenterQuadrature, ...
    'nCenterElements', projection.nCenterElements);
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

function [Mcenter, McenterInverse] = exactCenterMass(p)
%EXACTCENTERMASS Physical tensor-product DG mass and its block inverse.
% X is the faster center-coordinate index, hence M_X is the inner Kronecker
% factor. Element Jacobians are included because the adjoint is global.
massX = kron(spdiags(p.dg.X.jacobian(:), 0, ...
    p.dg.X.nElements, p.dg.X.nElements), p.dg.X.M);
massY = kron(spdiags(p.dg.Y.jacobian(:), 0, ...
    p.dg.Y.nElements, p.dg.Y.nElements), p.dg.Y.M);

inverseLocalX = sparse(inv(full(p.dg.X.M)));
inverseLocalY = sparse(inv(full(p.dg.Y.M)));
inverseMassX = kron(spdiags(1./p.dg.X.jacobian(:), 0, ...
    p.dg.X.nElements, p.dg.X.nElements), inverseLocalX);
inverseMassY = kron(spdiags(1./p.dg.Y.jacobian(:), 0, ...
    p.dg.Y.nElements, p.dg.Y.nElements), inverseLocalY);

Mcenter = sparse(kron(massY, massX));
McenterInverse = sparse(kron(inverseMassY, inverseMassX));
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

function projection = massProjectionData(p, params)
%MASSPROJECTIONDATA Quadrature-exact local multiplication in the DG basis.
%
% For nodal coefficients u_e on one rectangular element and a coefficient
% b evaluated at tensor-product Gauss-Lobatto points, the strong DG
% multiplication is
%
%   B_e(b) u_e = M_e^(-1) E' W diag(b) E u_e.
%
% E interpolates the local DG polynomial to quadrature points. The physical
% element Jacobian occurs in both M_e and W and cancels in the strong
% multiplication; it remains present in exactCenterMass for the adjoint.
nLocalX = p.dg.X.nLocal;
nLocalY = p.dg.Y.nLocal;
defaultQX = max(4, 2*nLocalX-1);
defaultQY = max(4, 2*nLocalY-1);
commonOrder = readParam(params, 'full2D_massQuadratureOrder', ...
    readParam(params, 'full2D_potentialQuadratureOrder', []));
if isempty(commonOrder)
    nQuadX = readParam(params, 'full2D_massQuadratureX', ...
        readParam(params, 'full2D_potentialQuadratureX', defaultQX));
    nQuadY = readParam(params, 'full2D_massQuadratureY', ...
        readParam(params, 'full2D_potentialQuadratureY', defaultQY));
else
    nQuadX = commonOrder;
    nQuadY = commonOrder;
end
nQuadX = max(nLocalX, round(nQuadX));
nQuadY = max(nLocalY, round(nQuadY));

[quadX, weightX] = JacobiGL(0, 0, nQuadX-1);
[quadY, weightY] = JacobiGL(0, 0, nQuadY-1);
interpolationX = Vandermonde1D(nLocalX-1, quadX)/p.dg.X.V;
interpolationY = Vandermonde1D(nLocalY-1, quadY)/p.dg.Y.V;
interpolation = kron(interpolationY, interpolationX);
weights = kron(weightY, weightX);
inverseMass = full(p.dg.MXY)\eye(size(p.dg.MXY));
centerProjection = (interpolation.*weights)*inverseMass.';

nQuadrature = numel(weights);
nLocal = nLocalX*nLocalY;
nElements = p.dg.X.nElements*p.dg.Y.nElements;
elementCenterIds = zeros(nLocal, nElements);
centerElementId = zeros(p.dg.nCenterDof, 1);
centerLocalId = zeros(p.dg.nCenterDof, 1);
centerQuadratureX = zeros(nQuadrature, nElements);
centerQuadratureY = zeros(nQuadrature, nElements);

element = 0;
for ey = 1:p.dg.Y.nElements
    yq = 0.5*(p.dg.Y.vertices(ey)+p.dg.Y.vertices(ey+1)) ...
        + p.dg.Y.jacobian(ey)*quadY(:);
    yIds = (ey-1)*nLocalY + (1:nLocalY);
    for ex = 1:p.dg.X.nElements
        element = element + 1;
        xq = 0.5*(p.dg.X.vertices(ex)+p.dg.X.vertices(ex+1)) ...
            + p.dg.X.jacobian(ex)*quadX(:);
        xIds = (ex-1)*nLocalX + (1:nLocalX);
        [localX, localY] = ndgrid(xIds, yIds);
        ids = sub2ind([p.dg.X.nDof, p.dg.Y.nDof], ...
            localX(:), localY(:));
        elementCenterIds(:, element) = ids;
        centerElementId(ids) = element;
        centerLocalId(ids) = (1:nLocal).';

        [Xq, Yq] = ndgrid(xq, yq);
        centerQuadratureX(:, element) = Xq(:);
        centerQuadratureY(:, element) = Yq(:);
    end
end

projection.nLocalCenter = nLocal;
projection.nCenterElements = nElements;
projection.nCenterQuadrature = nQuadrature;
projection.centerInterpolation = interpolation;
projection.centerProjection = centerProjection;
projection.elementCenterIds = elementCenterIds;
projection.centerElementId = centerElementId;
projection.centerLocalId = centerLocalId;
projection.centerQuadratureX = centerQuadratureX;
projection.centerQuadratureY = centerQuadratureY;
projection.quadratureOrder = [nQuadX, nQuadY];
targetEntries = max(1, floor(readParam(params, ...
    'full2D_variableMassApplyTargetEntries', 2e6)));
projection.elementChunkSize = min(nElements, max(1, ...
    floor(targetEntries/(p.relative.nDof*nQuadrature))));
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
        endpointInverseMassCorrection(mat, p, params, projection, ...
        axisName, referenceMass)
[inverseMass, gridX, gridY, fieldInfo] = inverseRelativeMassField( ...
    mat, p, axisName, referenceMass);
massInterpolant = griddedInterpolant({gridX, gridY}, inverseMass, ...
    'linear', 'nearest');

relativeX = repmat(p.relative.rhoX.cells(:), p.relative.NrhoY, 1);
relativeY = kron(p.relative.rhoY.cells(:), ...
    ones(p.relative.NrhoX, 1));

nRelative = p.relative.nDof;
nSamples = projection.nCenterQuadrature*projection.nCenterElements;
centerX = projection.centerQuadratureX(:);
centerY = projection.centerQuadratureY(:);
plusCorrection = zeros(nRelative, nSamples);
minusCorrection = zeros(nRelative, nSamples);
targetEntries = max(1, floor(readParam(params, ...
    'full2D_variableMassChunkTargetEntries', 2e6)));
sampleChunk = max(1, floor(targetEntries/nRelative));

for first = 1:sampleChunk:nSamples
    last = min(nSamples, first + sampleChunk - 1);
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

info = fieldInfo;
info.referenceRelativeMass = referenceMass;
info.inverseReferenceMass = 1/referenceMass;
info.plusCorrectionRange = rangeOf(plusCorrection);
info.minusCorrectionRange = rangeOf(minusCorrection);
relativeIds = reshape(1:p.relative.nDof, ...
    p.relative.NrhoX, p.relative.NrhoY);
reflectedIds = relativeIds(end:-1:1, end:-1:1);
reflectedMinus = minusCorrection(reflectedIds(:), :);
endpointScale = max(1, max(abs([plusCorrection(:); ...
    minusCorrection(:); 1/referenceMass])));
info.endpointExchangeRelativeDefect = max(abs(plusCorrection(:) ...
    - reflectedMinus(:)))/endpointScale;
info.endpointExchangeDefinition = ...
    'b_plus(R,rho) must equal b_minus(R,-rho).';
info.maxOuterXFaceCorrection = outerXFaceCorrection( ...
    massInterpolant, p, relativeX, relativeY, referenceMass);
info.outerXFacesMatchReference = info.maxOuterXFaceCorrection ...
    <= absoluteTolerance;
info.interpolation = ...
    'linear inverse mass at center Gauss-Lobatto points, nearest extension';
info.quadratureOrder = projection.quadratureOrder;
info.nProjectedSamples = nSamples;
info.chunkQuadratureSamples = sampleChunk;
end

function maximum = outerXFaceCorrection(massInterpolant, p, ...
        relativeX, relativeY, referenceMass)
%OUTERXFACECORRECTION Diagnose the homogeneous-lead assumption exactly on
%the two physical X faces. The projected volume samples do not generally
%coincide with every Y face node, so this check is deliberately separate.
faceY = p.dg.Y.nodes(:).';
faceX = [p.domain.X(1), p.domain.X(2)];
maximum = 0;
for side = 1:2
    X = faceX(side)*ones(size(faceY));
    plus = massInterpolant(X + 0.5*relativeX, ...
        faceY + 0.5*relativeY) - 1/referenceMass;
    minus = massInterpolant(X - 0.5*relativeX, ...
        faceY - 0.5*relativeY) - 1/referenceMass;
    maximum = max(maximum, max(abs([plus(:); minus(:)])));
end
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
        data.centerAdjointX, data.relativeDerivativeX, ...
        data.relativeAdjointX, data.massXPlus, data.massXMinus, ...
        data.coefficientX, data);
end
if data.hasYCorrection
    Y = Y + applyAxisCorrection(U, data.centerDerivativeY, ...
        data.centerAdjointY, data.relativeDerivativeY, ...
        data.relativeAdjointY, data.massYPlus, data.massYMinus, ...
        data.coefficientY, data);
end
y = Y(:);
end

function Y = applyAxisCorrection(U, C, CAdjoint, R, RAdjoint, ...
        massPlus, massMinus, coefficient, data)
plusGradient = 0.5*applyCenter(C, U) + R*U;
minusGradient = 0.5*applyCenter(C, U) - R*U;
plusFlux = applyMassProjection(data, plusGradient, massPlus, false);
minusFlux = applyMassProjection(data, minusGradient, massMinus, false);
plusAdjoint = 0.5*applyCenter(CAdjoint, plusFlux) ...
    + RAdjoint*plusFlux;
minusAdjoint = 0.5*applyCenter(CAdjoint, minusFlux) ...
    - RAdjoint*minusFlux;
Y = coefficient*(-plusAdjoint + minusAdjoint);
end

function Y = applyMassProjection(data, U, massSamples, absoluteBound)
%APPLYMASSPROJECTION Apply the local strong DG coefficient operator.
% massSamples(r,q,e) contains b(R_q,r)-b_ref. For absoluteBound=true the
% triangle-inequality operator |M^(-1)E'W| |b| |E| is applied; this is used
% only to construct a guaranteed row-sum bound for the rowabs preconditioner.
Y = zeros(size(U), 'like', U);
if absoluteBound
    interpolation = abs(data.centerInterpolation);
    projection = abs(data.centerProjection);
else
    interpolation = data.centerInterpolation;
    projection = data.centerProjection;
end

for first = 1:data.elementChunkSize:data.nCenterElements
    last = min(data.nCenterElements, first + data.elementChunkSize - 1);
    elements = first:last;
    centerIds = data.elementCenterIds(:, elements);
    nChunk = numel(elements);
    sampleIds = (first-1)*data.nCenterQuadrature+1: ...
        last*data.nCenterQuadrature;

    Ulocal = reshape(U(:, centerIds(:)), ...
        data.nRelative, data.nLocalCenter, nChunk);
    Uquadrature = pagemtimes(Ulocal, interpolation.');
    coefficients = reshape(massSamples(:, sampleIds), ...
        data.nRelative, data.nCenterQuadrature, nChunk);
    if absoluteBound
        coefficients = abs(coefficients);
    end
    localResidual = pagemtimes(coefficients.*Uquadrature, projection);
    Y(:, centerIds(:)) = reshape(localResidual, data.nRelative, []);
end
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
gpuData.centerAdjointX = gpuArray(data.centerAdjointX);
gpuData.centerAdjointY = gpuArray(data.centerAdjointY);
gpuData.relativeDerivativeX = gpuArray(data.relativeDerivativeX);
gpuData.relativeDerivativeY = gpuArray(data.relativeDerivativeY);
gpuData.relativeAdjointX = gpuArray(data.relativeAdjointX);
gpuData.relativeAdjointY = gpuArray(data.relativeAdjointY);
gpuData.centerInterpolation = gpuArray(data.centerInterpolation);
gpuData.centerProjection = gpuArray(data.centerProjection);
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
    correction = correction + assembleAxisCorrection(data, ...
        data.centerDerivativeX, data.centerAdjointX, ...
        data.relativeDerivativeX, data.relativeAdjointX, ...
        data.massXPlus, data.massXMinus, data.coefficientX);
end
if data.hasYCorrection
    correction = correction + assembleAxisCorrection(data, ...
        data.centerDerivativeY, data.centerAdjointY, ...
        data.relativeDerivativeY, data.relativeAdjointY, ...
        data.massYPlus, data.massYMinus, data.coefficientY);
end
correction = sparse(correction);
end

function A = assembleAxisCorrection(data, C, CAdjoint, R, RAdjoint, ...
        massPlus, massMinus, coefficient)
Icenter = speye(data.nCenter);
Irelative = speye(data.nRelative);
centerPart = 0.5*kron(C, Irelative);
relativePart = kron(Icenter, R);
Dplus = centerPart + relativePart;
Dminus = centerPart - relativePart;
centerAdjointPart = 0.5*kron(CAdjoint, Irelative);
relativeAdjointPart = kron(Icenter, RAdjoint);
DplusAdjoint = centerAdjointPart + relativeAdjointPart;
DminusAdjoint = centerAdjointPart - relativeAdjointPart;
Bplus = assembleMassMultiplication(data, massPlus);
Bminus = assembleMassMultiplication(data, massMinus);
A = coefficient*(-DplusAdjoint*Bplus*Dplus ...
    + DminusAdjoint*Bminus*Dminus);
A = sparse(A);
end

function B = assembleMassMultiplication(data, massSamples)
%ASSEMBLEMassMULTIPLICATION Explicit block-local projection for diagnostics.
% This path is only used when the complete transport matrix is requested.
rows = cell(data.nCenterElements, 1);
cols = cell(data.nCenterElements, 1);
values = cell(data.nCenterElements, 1);

for element = 1:data.nCenterElements
    localMatrix = spalloc(data.nRelative*data.nLocalCenter, ...
        data.nRelative*data.nLocalCenter, 0);
    firstSample = (element-1)*data.nCenterQuadrature;
    for iq = 1:data.nCenterQuadrature
        sample = firstSample + iq;
        relativeCoefficient = spdiags(massSamples(:, sample), 0, ...
            data.nRelative, data.nRelative);
        centerCoefficient = data.centerProjection(iq, :).' ...
            * data.centerInterpolation(iq, :);
        localMatrix = localMatrix ...
            + kron(sparse(centerCoefficient), relativeCoefficient);
    end

    centerIds = data.elementCenterIds(:, element);
    globalIds = (1:data.nRelative).' ...
        + (centerIds(:).'-1)*data.nRelative;
    globalIds = globalIds(:);
    [localRows, localCols, localValues] = find(localMatrix);
    rows{element} = globalIds(localRows);
    cols{element} = globalIds(localCols);
    values{element} = localValues;
end

B = sparse(vertcat(rows{:}), vertcat(cols{:}), vertcat(values{:}), ...
    data.nTotal, data.nTotal);
end

function diagonal = getCorrectionDiagonal(data)
diagonalMatrix = complex(zeros(data.nRelative, data.nCenter));
if data.hasXCorrection
    diagonalMatrix = diagonalMatrix + axisCorrectionDiagonal(data, ...
        data.centerDerivativeX, data.centerAdjointX, ...
        data.relativeDerivativeX, data.relativeAdjointX, ...
        data.massXPlus, data.massXMinus, data.coefficientX);
end
if data.hasYCorrection
    diagonalMatrix = diagonalMatrix + axisCorrectionDiagonal(data, ...
        data.centerDerivativeY, data.centerAdjointY, ...
        data.relativeDerivativeY, data.relativeAdjointY, ...
        data.massYPlus, data.massYMinus, data.coefficientY);
end
diagonal = diagonalMatrix(:);
end

function diagonal = axisCorrectionDiagonal(data, C, CAdjoint, R, ...
        RAdjoint, plusMass, minusMass, coefficient)
diagonal = complex(zeros(data.nRelative, data.nCenter));
relativeKernel = RAdjoint.*R.';
diagonalR = full(diag(R));
diagonalRAdjoint = full(diag(RAdjoint));
for centerId = 1:data.nCenter
    plusDiagonal = derivativeBlockDiagonal(data, C, CAdjoint, ...
        relativeKernel, diagonalR, diagonalRAdjoint, plusMass, ...
        centerId, 1);
    minusDiagonal = derivativeBlockDiagonal(data, C, CAdjoint, ...
        relativeKernel, diagonalR, diagonalRAdjoint, minusMass, ...
        centerId, -1);
    diagonal(:, centerId) = coefficient*(-plusDiagonal + minusDiagonal);
end
end

function diagonal = derivativeBlockDiagonal(data, C, CAdjoint, ...
        relativeKernel, diagonalR, diagonalRAdjoint, massSamples, ...
        centerId, signR)
[centerDiagonal, leftCross, rightCross, centerRoundTrip] = ...
    centerBlockCoefficients(data, C, CAdjoint, massSamples, centerId);
diagonal = relativeKernel*centerDiagonal ...
    + 0.5*signR*(leftCross.*diagonalR ...
    + diagonalRAdjoint.*rightCross) ...
    + 0.25*centerRoundTrip;
end

function rowAbs = getCorrectionRowAbsBound(data)
rowAbsMatrix = zeros(data.nRelative, data.nCenter);
if data.hasXCorrection
    rowAbsMatrix = rowAbsMatrix + axisCorrectionRowAbsBound(data, ...
        data.centerDerivativeX, data.centerAdjointX, ...
        data.relativeDerivativeX, data.relativeAdjointX, ...
        data.massXPlus, data.massXMinus, data.coefficientX);
end
if data.hasYCorrection
    rowAbsMatrix = rowAbsMatrix + axisCorrectionRowAbsBound(data, ...
        data.centerDerivativeY, data.centerAdjointY, ...
        data.relativeDerivativeY, data.relativeAdjointY, ...
        data.massYPlus, data.massYMinus, data.coefficientY);
end
rowAbs = rowAbsMatrix(:);
end

function bound = axisCorrectionRowAbsBound(data, C, CAdjoint, R, ...
        RAdjoint, plusMass, minusMass, coefficient)
absC = abs(C);
absR = abs(R);
rowMagnitude = repmat(full(sum(absR, 2)), 1, data.nCenter) ...
    + 0.5*repmat(full(sum(absC, 2)).', data.nRelative, 1);
plusWeighted = applyMassProjection(data, rowMagnitude, plusMass, true);
minusWeighted = applyMassProjection(data, rowMagnitude, minusMass, true);
plusBound = abs(RAdjoint)*plusWeighted ...
    + 0.5*applyCenter(abs(CAdjoint), plusWeighted);
minusBound = abs(RAdjoint)*minusWeighted ...
    + 0.5*applyCenter(abs(CAdjoint), minusWeighted);
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
    block = block + axisCorrectionBlock(data, data.centerDerivativeX, ...
        data.centerAdjointX, data.relativeDerivativeX, ...
        data.relativeAdjointX, data.massXPlus, data.massXMinus, ...
        data.coefficientX, centerId);
end
if data.hasYCorrection
    block = block + axisCorrectionBlock(data, data.centerDerivativeY, ...
        data.centerAdjointY, data.relativeDerivativeY, ...
        data.relativeAdjointY, data.massYPlus, data.massYMinus, ...
        data.coefficientY, centerId);
end
block = sparse(block);
end

function block = axisCorrectionBlock(data, C, CAdjoint, R, RAdjoint, ...
        plusMass, minusMass, coefficient, centerId)
plusBlock = derivativeCenterBlock(data, C, CAdjoint, R, RAdjoint, ...
    plusMass, centerId, 1);
minusBlock = derivativeCenterBlock(data, C, CAdjoint, R, RAdjoint, ...
    minusMass, centerId, -1);
block = coefficient*(-plusBlock + minusBlock);
end

function block = derivativeCenterBlock(data, C, CAdjoint, R, RAdjoint, ...
        massSamples, centerId, signR)
[centerDiagonal, leftCross, rightCross, centerRoundTrip] = ...
    centerBlockCoefficients(data, C, CAdjoint, massSamples, centerId);
block = RAdjoint*spdiags(centerDiagonal, 0, data.nRelative, ...
    data.nRelative)*R ...
    + 0.5*signR*(spdiags(leftCross, 0, data.nRelative, ...
    data.nRelative)*R + RAdjoint*spdiags(rightCross, 0, ...
    data.nRelative, data.nRelative)) ...
    + 0.25*spdiags(centerRoundTrip, 0, data.nRelative, ...
    data.nRelative);
block = sparse(block);
end

function [centerDiagonal, leftCross, rightCross, centerRoundTrip] = ...
        centerBlockCoefficients(data, C, CAdjoint, massSamples, centerId)
%CENTERBLOCKCOEFFICIENTS Exact rho-block coefficients of B_e(b).
% Although B_e is dense inside one DG element, it is diagonal in the
% relative-coordinate index. These four projected coefficient vectors are
% sufficient to form the exact center-coordinate diagonal block used by
% block Jacobi without assembling the complete four-dimensional matrix.
element = data.centerElementId(centerId);
localId = data.centerLocalId(centerId);
centerIds = data.elementCenterIds(:, element);
sampleIds = (element-1)*data.nCenterQuadrature ...
    + (1:data.nCenterQuadrature);
massElement = massSamples(:, sampleIds);

leftAtQuadrature = data.centerProjection ...
    * full(CAdjoint(centerId, centerIds)).';
rightAtQuadrature = data.centerInterpolation ...
    * full(C(centerIds, centerId));
trialAtCenter = data.centerInterpolation(:, localId);
testAtCenter = data.centerProjection(:, localId);

weights = [testAtCenter.*trialAtCenter, ...
    leftAtQuadrature.*trialAtCenter, ...
    testAtCenter.*rightAtQuadrature];
coefficients = massElement*weights;
centerDiagonal = coefficients(:, 1);
leftCross = coefficients(:, 2);
rightCross = coefficients(:, 3);

% The center-center path C^* B C is different from the two cross terms:
% C may first send the trial value into a neighbouring DG element, B acts
% locally in that neighbour, and C^* sends the result back. Therefore all
% elements touched by both the C column and the C^* row must contribute;
% restricting this term to the element that owns centerId would omit face
% couplings from the Jacobi diagonal and the rho-block preconditioner.
columnSupport = find(C(:, centerId));
rowSupport = find(CAdjoint(centerId, :)).';
columnElements = unique(data.centerElementId(columnSupport));
rowElements = unique(data.centerElementId(rowSupport));
roundTripElements = intersect(columnElements, rowElements);
centerRoundTrip = zeros(data.nRelative, 1);
for id = 1:numel(roundTripElements)
    roundTripElement = roundTripElements(id);
    roundTripCenterIds = data.elementCenterIds(:, roundTripElement);
    roundTripSampleIds = (roundTripElement-1)*data.nCenterQuadrature ...
        + (1:data.nCenterQuadrature);
    left = data.centerProjection ...
        * full(CAdjoint(centerId, roundTripCenterIds)).';
    right = data.centerInterpolation ...
        * full(C(roundTripCenterIds, centerId));
    centerRoundTrip = centerRoundTrip ...
        + massSamples(:, roundTripSampleIds)*(left.*right);
end
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
