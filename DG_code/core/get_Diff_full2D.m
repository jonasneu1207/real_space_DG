function [A, rhs, info] = get_Diff_full2D(mat, p, boundary)
%GET_DIFF_FULL2D Full-2D DG/FV transport operator in rho basis.
%
% This is the first assembled part of the full physical 2D Wigner equation.
% The unknown is stored in the relative-coordinate basis
%
%   F = F(X, Y, rho_x, rho_y).
%
% Global DOF order:
%   F(iX, iY, iRhoX, iRhoY) is vectorized as F(:). Hence X-DG DOFs are the
%   fastest index, followed by Y-DG DOFs, rho_x FV cells and rho_y FV cells.
%   If C acts on the center coordinates (X,Y) and R acts on the relative
%   coordinates (rho_x,rho_y), the lifted matrix is
%
%       kron(R, C).
%
% Transport model used here:
%   A_X = -1i*d/d(rho_x),  A_Y = -1i*d/d(rho_y)
%
% and the center-coordinate DG transport operator is assembled as
%
%   Qx * (A_X*d/dX) + Qy * (A_Y*d/dY).
%
% The volume and central numerical flux terms stay sparse because A_X and
% A_Y are finite-volume derivative matrices in rho. Interior element faces
% use the scalar Rusanov/Lax-Friedrichs approximation
%
%   |A_n| approx alphaLF * I,
%
% so the stabilizing face term also remains sparsity-friendly. Dense
% characteristic splits are used only on physical boundary faces to impose
% Source/Drain inflow and specular reflection cleanly.

if nargin == 1
    p = mat;
    mat = struct;
    boundary = [];
elseif nargin < 3
    boundary = [];
end

if isempty(boundary)
    if isempty(fieldnames(mat))
        error('DG:Full2D:MissingBoundaryData', ...
            'get_Diff_full2D needs mat and boundary data for Source/Drain RHS assembly.');
    end
    boundary = get_Boundary_full2D(mat, p, p.EfL, p.EfR, []);
end

params = getDGParams(mat);
fluxType = lower(readParam(params, 'full2D_fluxType', ...
    readParam(params, 'fluxType', 'rusanov')));
if any(strcmpi(fluxType, {'lf', 'lax-friedrichs', 'rusanov-local'}))
    fluxType = 'rusanov';
end

axisOpsX = oneDimensionalDGOperators(p.dg.X);
axisOpsY = oneDimensionalDGOperators(p.dg.Y);
centerOps = liftCenterOperators(p, axisOpsX, axisOpsY);
relativeOps = relativeTransportOperators(p, params, fluxType);
[scaleX, scaleY, scaleInfo] = transportScales(mat, p, params);

operatorParts = {};
rhsParts = {};

% X transport: central DG volume/face contribution plus optional Rusanov
% stabilization on interior vertical element interfaces.
operatorParts = addOperatorPart(operatorParts, 'X central rho transport', ...
    scaleX, centerOps.X.H1, p.relative.Ax);
operatorParts = addInteriorRusanovPart(operatorParts, 'X interior Rusanov', ...
    scaleX, centerOps.X.H2Interior, relativeOps.X.AabsLF, relativeOps.X.betaLF);

% Characteristic boundary penalty on the outer X faces. This is the only
% place in the X direction where dense characteristic matrices enter.
operatorParts = addOperatorPart(operatorParts, 'X physical boundary characteristic penalty', ...
    scaleX, centerOps.X.H2Boundary, relativeOps.X.AabsChar);

% Masked Source/Drain face pieces are not reservoirs. The base boundary
% penalty enforces incoming characteristics against rhoBoundary; on closed
% pieces we replace this by the specular condition rho_in = R_x*rho_inside.
[operatorParts, rhsParts] = addSourceBoundary(operatorParts, rhsParts, ...
    boundary.physical.XLeft, centerOps.X.leftLiftByY, ...
    centerOps.X.leftVectorByY, relativeOps.X, scaleX);
[operatorParts, rhsParts] = addDrainBoundary(operatorParts, rhsParts, ...
    boundary.physical.XRight, centerOps.X.rightLiftByY, ...
    centerOps.X.rightVectorByY, relativeOps.X, scaleX);

% Y transport: same structure as X, but there are no reservoirs on the
% bottom/top physical boundaries. The incoming characteristic data there are
% reflected with R_y.
operatorParts = addOperatorPart(operatorParts, 'Y central rho transport', ...
    scaleY, centerOps.Y.H1, p.relative.Ay);
operatorParts = addInteriorRusanovPart(operatorParts, 'Y interior Rusanov', ...
    scaleY, centerOps.Y.H2Interior, relativeOps.Y.AabsLF, relativeOps.Y.betaLF);
operatorParts = addOperatorPart(operatorParts, 'Y physical boundary characteristic penalty', ...
    scaleY, centerOps.Y.H2Boundary, relativeOps.Y.AabsChar);
operatorParts = addYReflection(operatorParts, boundary.physical.YBottom, ...
    centerOps.Y.bottomLift, relativeOps.Y.AplusChar, -scaleY);
operatorParts = addYReflection(operatorParts, boundary.physical.YTop, ...
    centerOps.Y.topLift, relativeOps.Y.AminusChar, scaleY);

assembleMatrix = shouldAssembleMatrix(params, p.index.nTotal);
if assembleMatrix
    A = assembleKronOperator(operatorParts, p.dg.nCenterDof, p.relative.nDof);
    rhs = assembleRhs(rhsParts, p.dg.nCenterDof, p.relative.nDof);
else
    A = [];
    rhs = [];
end

info = struct;
info.full2D = true;
info.fullMatrixAssembled = assembleMatrix;
info.matrixFreeAvailable = true;
info.reason = matrixReason(assembleMatrix, params, p.index.nTotal);
info.fluxType = fluxType;
info.dofOrder = p.index.order;
info.totalDof = p.index.nTotal;
info.centerDof = p.dg.nCenterDof;
info.relativeDof = p.relative.nDof;
info.transportScales = scaleInfo;
info.alphaLF = [relativeOps.X.alphaLF, relativeOps.Y.alphaLF];
info.thetaLF = [relativeOps.X.thetaLF, relativeOps.Y.thetaLF];
info.betaLF = [relativeOps.X.betaLF, relativeOps.Y.betaLF];
info.boundary.sourceContactDof = boundary.physical.XLeft.nContactFaceDof;
info.boundary.sourceClosedDof = boundary.physical.XLeft.nClosedFaceDof;
info.boundary.drainContactDof = boundary.physical.XRight.nContactFaceDof;
info.boundary.drainClosedDof = boundary.physical.XRight.nClosedFaceDof;
info.operatorPartSummary = summarizeOperatorParts(operatorParts);
info.rhsPartSummary = summarizeRhsParts(rhsParts);
info.apply = @(u) applyKronOperator(operatorParts, u, ...
    p.dg.nCenterDof, p.relative.nDof);
info.assemble = @() assembleKronOperator(operatorParts, ...
    p.dg.nCenterDof, p.relative.nDof);
info.assembleRhs = @() assembleRhs(rhsParts, ...
    p.dg.nCenterDof, p.relative.nDof);
info.getDiagonal = @() assembleKronDiagonal(operatorParts, ...
    p.dg.nCenterDof, p.relative.nDof);
end

function axisOps = oneDimensionalDGOperators(axis)
%ONEDIMENSIONALDGOPERATORS Build 1D DG matrices for one rectangular axis.
%
% The formulas mirror the established 1D rho solver. H1 contains the
% derivative and the central numerical flux contribution. H2 multiplies an
% absolute-value flux matrix. We split H2 into an interior part, where the
% sparse scalar Rusanov approximation is used, and a physical-boundary part,
% where characteristic matrices are allowed.

nLocal = axis.nLocal;
nElement = axis.nElements;
nDof = axis.nDof;
invM = sparse(inv(full(axis.M)));
stiffness = sparse(axis.S);
jacobian = axis.jacobian(:);

P1 = sparse(1, 1, 1, nLocal, nLocal);
P2 = sparse(1, nLocal, 1, nLocal, nLocal);
P3 = sparse(nLocal, nLocal, 1, nLocal, nLocal);
P4 = sparse(nLocal, 1, 1, nLocal, nLocal);

K1 = invM*(stiffness + 0.5*P1 - 0.5*P3);
K2 = invM*(0.5*P1 + 0.5*P3);
K2L = invM*(0.5*P1);
K2R = invM*(0.5*P3);
K3 = -0.5*invM*P2;
K4 = K3;
K5 = 0.5*invM*P4;
K6 = -K5;

upperScale = spdiags(1./jacobian, 1, nElement, nElement);
lowerScale = spdiags(1./jacobian(2:end), -1, nElement, nElement);
diagScale = spdiags(1./jacobian, 0, nElement, nElement);

H1 = kron(upperScale, K5) + kron(lowerScale, K3) + kron(diagScale, K1);
H2 = kron(upperScale, K6) + kron(lowerScale, K4) + kron(diagScale, K2);

leftElement = sparse(1, 1, 1, nElement, nElement);
rightElement = sparse(nElement, nElement, 1, nElement, nElement);
H2LeftBoundary = kron(leftElement*spdiags(1./jacobian, 0, nElement, nElement), K2L);
H2RightBoundary = kron(rightElement*spdiags(1./jacobian, 0, nElement, nElement), K2R);
H2Boundary = H2LeftBoundary + H2RightBoundary;

leftTrace = sparse(1, 1, 1, nDof, 1);
rightTrace = sparse(nDof, 1, 1, nDof, 1);
leftLift = sparse(nDof, 1);
rightLift = sparse(nDof, 1);
leftLift(1:nLocal) = invM(:, 1)/jacobian(1);
rightLift(nDof-nLocal+1:nDof) = invM(:, nLocal)/jacobian(end);

axisOps = struct;
axisOps.H1 = sparse(H1);
axisOps.H2 = sparse(H2);
axisOps.H2Boundary = sparse(H2Boundary);
axisOps.H2Interior = sparse(H2 - H2Boundary);
axisOps.leftBoundary = sparse(leftLift*leftTrace.');
axisOps.rightBoundary = sparse(rightLift*rightTrace.');
axisOps.leftLift = leftLift;
axisOps.rightLift = rightLift;
axisOps.nDof = nDof;
end

function centerOps = liftCenterOperators(p, axisOpsX, axisOpsY)
%LIFTCENTEROPERATORS Lift 1D DG operators to the rectangular X-Y grid.
%
% Center-coordinate ordering is (X,Y), with X fastest. Thus an X-only
% operator is kron(I_Y, O_X), while a Y-only operator is kron(O_Y, I_X).

IX = speye(p.dg.X.nDof);
IY = speye(p.dg.Y.nDof);

centerOps = struct;
centerOps.X.H1 = kron(IY, axisOpsX.H1);
centerOps.X.H2Interior = kron(IY, axisOpsX.H2Interior);
centerOps.X.H2Boundary = kron(IY, axisOpsX.H2Boundary);
centerOps.X.leftLiftByY = @(iyMask) kron(spdiags(double(iyMask(:)), 0, ...
    p.dg.Y.nDof, p.dg.Y.nDof), axisOpsX.leftBoundary);
centerOps.X.rightLiftByY = @(iyMask) kron(spdiags(double(iyMask(:)), 0, ...
    p.dg.Y.nDof, p.dg.Y.nDof), axisOpsX.rightBoundary);
centerOps.X.leftVectorByY = @(iy) kron(unitVector(p.dg.Y.nDof, iy), axisOpsX.leftLift);
centerOps.X.rightVectorByY = @(iy) kron(unitVector(p.dg.Y.nDof, iy), axisOpsX.rightLift);

centerOps.Y.H1 = kron(axisOpsY.H1, IX);
centerOps.Y.H2Interior = kron(axisOpsY.H2Interior, IX);
centerOps.Y.H2Boundary = kron(axisOpsY.H2Boundary, IX);
centerOps.Y.bottomLift = kron(axisOpsY.leftBoundary, IX);
centerOps.Y.topLift = kron(axisOpsY.rightBoundary, IX);
end

function relativeOps = relativeTransportOperators(p, params, fluxType)
%RELATIVETRANSPORTOPERATORS Characteristic and Rusanov matrices in rho.
%
% A_X and A_Y themselves are sparse finite-volume derivative operators. The
% characteristic splits A^+/A^- are only used on physical boundary faces.
% Interior Rusanov stabilization uses betaLF*I instead of a dense |A|.

tol = readParam(params, 'full2D_characteristicTol', 1e-12);
IX = speye(p.relative.NrhoX);
IY = speye(p.relative.NrhoY);
IR = p.relative.identity;

[AplusX1D, AminusX1D, lambdaX] = splitCharacteristicMatrix(p.relative.Ax1D, tol);
[AplusY1D, AminusY1D, lambdaY] = splitCharacteristicMatrix(p.relative.Ay1D, tol);

thetaDefault = 0.01;
if strcmpi(fluxType, 'central')
    thetaX = 0;
    thetaY = 0;
elseif strcmpi(fluxType, 'rusanov')
    thetaX = readParam(params, 'full2D_thetaLF_X', ...
        readParam(params, 'full2D_thetaLF', readParam(params, 'thetaLF', thetaDefault)));
    thetaY = readParam(params, 'full2D_thetaLF_Y', ...
        readParam(params, 'full2D_thetaLF', readParam(params, 'thetaLF', thetaDefault)));
else
    error('DG:Full2D:UnknownFluxType', ...
        'Unknown full2D_fluxType "%s". Use "rusanov" or "central".', fluxType);
end

alphaX = readParam(params, 'full2D_alphaLF_X', max(abs(lambdaX)));
alphaY = readParam(params, 'full2D_alphaLF_Y', max(abs(lambdaY)));

relativeOps = struct;
relativeOps.X.AplusChar = sparse(kron(IY, AplusX1D));
relativeOps.X.AminusChar = sparse(kron(IY, AminusX1D));
relativeOps.X.AabsChar = sparse(relativeOps.X.AplusChar - relativeOps.X.AminusChar);
relativeOps.X.alphaLF = alphaX;
relativeOps.X.thetaLF = thetaX;
relativeOps.X.betaLF = thetaX*alphaX;
relativeOps.X.AabsLF = sparse(relativeOps.X.betaLF*IR);

relativeOps.Y.AplusChar = sparse(kron(AplusY1D, IX));
relativeOps.Y.AminusChar = sparse(kron(AminusY1D, IX));
relativeOps.Y.AabsChar = sparse(relativeOps.Y.AplusChar - relativeOps.Y.AminusChar);
relativeOps.Y.alphaLF = alphaY;
relativeOps.Y.thetaLF = thetaY;
relativeOps.Y.betaLF = thetaY*alphaY;
relativeOps.Y.AabsLF = sparse(relativeOps.Y.betaLF*IR);
end

function [operatorParts, rhsParts] = addSourceBoundary(operatorParts, rhsParts, ...
        source, leftLiftByY, leftVectorByY, relativeX, scaleX)
contactMask = boundaryContactMask(source);
closedMask = ~contactMask;
if any(closedMask)
    closedCenter = leftLiftByY(closedMask);
    reflectedInflow = relativeX.AplusChar*source.nonContactGhostOperator;
    operatorParts = addOperatorPart(operatorParts, ...
        'X-left closed-face reflection correction', ...
        -scaleX, closedCenter, reflectedInflow);
end

contactIds = find(contactMask);
for id = contactIds(:).'
    rhsParts = addRhsPart(rhsParts, sprintf('X-left Source inflow yDof=%d', id), ...
        scaleX, leftVectorByY(id), ...
        relativeX.AplusChar*source.rhoBoundary(:, id));
end
end

function [operatorParts, rhsParts] = addDrainBoundary(operatorParts, rhsParts, ...
        drain, rightLiftByY, rightVectorByY, relativeX, scaleX)
contactMask = boundaryContactMask(drain);
closedMask = ~contactMask;
if any(closedMask)
    closedCenter = rightLiftByY(closedMask);
    reflectedInflow = relativeX.AminusChar*drain.nonContactGhostOperator;
    operatorParts = addOperatorPart(operatorParts, ...
        'X-right closed-face reflection correction', ...
        scaleX, closedCenter, reflectedInflow);
end

contactIds = find(contactMask);
for id = contactIds(:).'
    rhsParts = addRhsPart(rhsParts, sprintf('X-right Drain inflow yDof=%d', id), ...
        scaleX, rightVectorByY(id), ...
        -relativeX.AminusChar*drain.rhoBoundary(:, id));
end
end

function operatorParts = addYReflection(operatorParts, side, centerLift, charMatrix, coefficient)
%ADDYREFLECTION Replace reservoir data on Y boundaries by reflected ghosts.
%
% The base physical-boundary penalty uses the incoming characteristic matrix
% times the interior trace. On a specular wall the incoming state is not an
% external reservoir; it is R_y times the interior state. Therefore we add
% the correction that changes rho_in = 0 into rho_in = R_y*rho_inside.

operatorParts = addOperatorPart(operatorParts, ...
    [side.face, ' specular reflection correction'], ...
    coefficient, centerLift, charMatrix*side.R_y);
end

function mask = boundaryContactMask(side)
if isfield(side, 'contactMask') && ~isempty(side.contactMask)
    mask = logical(side.contactMask(:));
else
    mask = true(size(side.rhoBoundary, 2), 1);
end
end

function operatorParts = addInteriorRusanovPart(operatorParts, name, ...
        coefficient, centerMatrix, relativeMatrix, betaLF)
if betaLF == 0
    return
end
operatorParts = addOperatorPart(operatorParts, name, ...
    coefficient, centerMatrix, relativeMatrix);
end

function operatorParts = addOperatorPart(operatorParts, name, coefficient, ...
        centerMatrix, relativeMatrix)
part = struct;
part.name = name;
part.coefficient = coefficient;
part.centerMatrix = sparse(centerMatrix);
part.relativeMatrix = sparse(relativeMatrix);
part.centerNnz = nnz(part.centerMatrix);
part.relativeNnz = nnz(part.relativeMatrix);
operatorParts{end+1} = part;
end

function rhsParts = addRhsPart(rhsParts, name, coefficient, centerVector, relativeVector)
part = struct;
part.name = name;
part.coefficient = coefficient;
part.centerVector = sparse(centerVector(:));
part.relativeVector = sparse(relativeVector(:));
part.centerNnz = nnz(part.centerVector);
part.relativeNnz = nnz(part.relativeVector);
rhsParts{end+1} = part;
end

function A = assembleKronOperator(operatorParts, nCenter, nRelative)
A = spalloc(nCenter*nRelative, nCenter*nRelative, 0);
for ip = 1:numel(operatorParts)
    part = operatorParts{ip};
    A = A + part.coefficient*kron(part.relativeMatrix, part.centerMatrix);
end
A = sparse(A);
end

function diagonal = assembleKronDiagonal(operatorParts, nCenter, nRelative)
%ASSEMBLEKRONDIAGONAL Diagonal of the matrix-free transport operator.
%
% For the global ordering (center fastest, relative slowest), each separable
% contribution coefficient*kron(R,C) has diagonal
%
%   coefficient*kron(diag(R), diag(C)).
%
% This is the quantity used by the Full-2D Jacobi preconditioner; no global
% sparse matrix has to be formed.

diagonal = complex(zeros(nCenter*nRelative, 1));
for ip = 1:numel(operatorParts)
    part = operatorParts{ip};
    centerDiagonal = full(diag(part.centerMatrix));
    relativeDiagonal = full(diag(part.relativeMatrix));
    diagonal = diagonal ...
        + part.coefficient*kron(relativeDiagonal, centerDiagonal);
end
end

function rhs = assembleRhs(rhsParts, nCenter, nRelative)
rhs = sparse(nCenter*nRelative, 1);
for ip = 1:numel(rhsParts)
    part = rhsParts{ip};
    rhs = rhs + part.coefficient*kron(part.relativeVector, part.centerVector);
end
rhs = full(rhs);
end

function y = applyKronOperator(operatorParts, u, nCenter, nRelative)
if numel(u) ~= nCenter*nRelative
    error('DG:Full2D:InvalidOperatorInput', ...
        'Input has %d entries, expected %d.', numel(u), nCenter*nRelative);
end
U = reshape(u, nCenter, nRelative);
Y = zeros(nCenter, nRelative);
for ip = 1:numel(operatorParts)
    part = operatorParts{ip};
    Y = Y + part.coefficient*(part.centerMatrix*U*part.relativeMatrix.');
end
y = Y(:);
end

function tf = shouldAssembleMatrix(params, nTotal)
mode = readParam(params, 'full2D_assembleDiffMatrix', 'auto');
maxDof = readParam(params, 'full2D_maxAssembledDof', 50000);
if islogical(mode)
    tf = mode;
elseif isnumeric(mode)
    tf = mode ~= 0;
else
    switch lower(char(mode))
        case {'auto'}
            tf = nTotal <= maxDof;
        case {'true', 'yes', 'on', 'assemble', 'always'}
            tf = true;
        case {'false', 'no', 'off', 'matrix-free', 'never'}
            tf = false;
        otherwise
            error('DG:Full2D:UnknownAssemblyMode', ...
                'Unknown full2D_assembleDiffMatrix mode "%s".', char(mode));
    end
end
end

function reason = matrixReason(assembled, params, nTotal)
mode = readParam(params, 'full2D_assembleDiffMatrix', 'auto');
if assembled
    reason = sprintf('Diff matrix assembled in %s mode for %d total DOFs.', ...
        char(string(mode)), nTotal);
else
    reason = sprintf(['Diff matrix kept matrix-free in %s mode for %d total DOFs. ', ...
        'Use info.apply(u), info.assemble(), or increase full2D_maxAssembledDof for small tests.'], ...
        char(string(mode)), nTotal);
end
end

function summary = summarizeOperatorParts(operatorParts)
summary = struct('name', {}, 'coefficient', {}, 'centerNnz', {}, 'relativeNnz', {});
for ip = 1:numel(operatorParts)
    summary(ip).name = operatorParts{ip}.name;
    summary(ip).coefficient = operatorParts{ip}.coefficient;
    summary(ip).centerNnz = operatorParts{ip}.centerNnz;
    summary(ip).relativeNnz = operatorParts{ip}.relativeNnz;
end
end

function summary = summarizeRhsParts(rhsParts)
summary = struct('name', {}, 'coefficient', {}, 'centerNnz', {}, 'relativeNnz', {});
for ip = 1:numel(rhsParts)
    summary(ip).name = rhsParts{ip}.name;
    summary(ip).coefficient = rhsParts{ip}.coefficient;
    summary(ip).centerNnz = rhsParts{ip}.centerNnz;
    summary(ip).relativeNnz = rhsParts{ip}.relativeNnz;
end
end

function [Aplus, Aminus, lambda] = splitCharacteristicMatrix(A, tol)
[Phi, Lambda] = eig(full(A));
lambda = diag(Lambda);
lambdaPlus = lambda;
lambdaMinus = lambda;
lambdaPlus(real(lambdaPlus) <= tol) = 0;
lambdaMinus(real(lambdaMinus) >= -tol) = 0;
PhiInv = Phi\eye(size(Phi, 1));
Aplus = Phi*diag(lambdaPlus)*PhiInv;
Aminus = Phi*diag(lambdaMinus)*PhiInv;
Aplus(abs(Aplus) < 1e-12) = 0;
Aminus(abs(Aminus) < 1e-12) = 0;
end

function [scaleX, scaleY, info] = transportScales(mat, p, params)
constants = physicalConstants;
defaultMeX = readField(mat, 'me_x_ch', 0.041);
defaultMeY = readField(mat, 'me_y_ch', defaultMeX);
defaultScaleX = constants.hbar/(massToKg(defaultMeX));
defaultScaleY = constants.hbar/(massToKg(defaultMeY));

scaleX = readParam(params, 'full2D_Q_diff_x', ...
    readParam(params, 'full2D_transportScaleX', defaultScaleX));
scaleY = readParam(params, 'full2D_Q_diff_y', ...
    readParam(params, 'full2D_transportScaleY', defaultScaleY));
scaleFactor = readParam(params, 'full2D_diff_scale', ...
    readParam(params, 'rho_diff_scale', 1));
scaleX = scaleFactor*scaleX;
scaleY = scaleFactor*scaleY;

info = struct;
info.Qx = scaleX;
info.Qy = scaleY;
info.scaleFactor = scaleFactor;
info.note = ['The current sparse Kron implementation uses scalar reference ', ...
    'transport scales. Spatially varying effective mass is a later extension.'];
info.dofOrder = p.index.order;
end

function massKg = massToKg(massValue)
constants = physicalConstants;
if abs(massValue) > 1e-25
    massKg = massValue*constants.m0;
else
    massKg = massValue;
end
if massKg <= 0
    massKg = 0.041*constants.m0;
end
end

function e = unitVector(n, id)
e = sparse(id, 1, 1, n, 1);
end

function params = getDGParams(mat)
params = struct;
if isfield(mat, 'dg') && isfield(mat.dg, 'params')
    params = mat.dg.params;
end
end

function value = readParam(params, name, defaultValue)
value = defaultValue;
if isstruct(params)
    if isfield(params, name) && ~isempty(params.(name))
        value = params.(name);
    elseif isfield(params, 'full2D') && isstruct(params.full2D) ...
            && isfield(params.full2D, name) && ~isempty(params.full2D.(name))
        value = params.full2D.(name);
    end
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
