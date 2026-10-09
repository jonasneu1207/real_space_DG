function [G, info] = get_Drift_full2D(mat, p, Vxy)
%GET_DRIFT_FULL2D Potential and CAP operator for Full-2D Wigner-DG.
%
% Unknown and vector order:
%   F(iRhoX,iRhoY,iX,iY) is stored as F(:). rho_x is the fastest
%   coordinate, followed by rho_y, the X-DG DOF and the Y-DG DOF.
%
% The default FV-consistent potential path evaluates
%
%   DeltaV = V(X+rho_x/2,Y+rho_y/2) - V(X-rho_x/2,Y-rho_y/2)
%
% at tensor-product FV vertices in (rho_x,rho_y). If P maps relative FV
% cell values to the surrounding vertices by summation, the mass-normalized
% relative-coordinate potential block is
%
%   C_V = (1/16) * P' * diag(DeltaV) * P.
%
% This is the direct two-dimensional analogue of the 1D rho/FV formula
% (1/4)*P'*diag(DeltaV)*P. It couples at most the 3-by-3 neighborhood of a
% relative FV cell and therefore remains sparse. C_V is projected over the
% tensor-product X/Y DG basis element by element.
%
% The legacy collocated path is retained for controlled comparisons. CAP is
% diagonal and separable in rho_x/rho_y in both paths. Physical boundary
% fluxes are deliberately absent here: Source/Drain characteristic splits,
% X-boundary eigendecompositions and all numerical fluxes are assembled
% exclusively by get_Diff_full2D.

if nargin < 3 || isempty(Vxy)
    if isfield(mat, 'V') && ~isempty(mat.V)
        Vxy = mat.V;
    else
        Vxy = zeros(mat.Nx, mat.Ny);
    end
end

params = getDGParams(mat);
potentialMode = normalizePotentialDiscretization(readParam(params, ...
    'full2D_potentialDiscretization', 'fv-consistent'));
cap = get_CAP_full2D(p);
[Vfield, gridX, gridY, potentialInfo] = potentialOnMaterialGrid(mat, p, Vxy);

data = driftData(mat, p, params, cap, Vfield, gridX, gridY, potentialMode);
assembleMatrix = shouldAssembleDriftMatrix(params, p.index.nTotal);

switch potentialMode
    case 'collocated'
        storeDiagonal = assembleMatrix ...
            || shouldStoreDriftDiagonal(params, p.index.nTotal);
        if storeDiagonal
            data.storedDiagonal = materializeCollocatedDiagonal(data);
        end
    case 'fv-consistent'
        data = addFVConsistentData(data, p, params);
        storePotential = assembleMatrix || shouldStorePotentialDifference( ...
            params, data.nPotentialSamples);
        if storePotential
            data.storedPotentialDifference = ...
                materializePotentialDifference(data);
        end
end

if assembleMatrix
    G = assembleDriftMatrix(data);
else
    G = [];
end

info = struct;
info.full2D = true;
info.fullPotentialAssembled = assembleMatrix;
info.matrixFreeAvailable = true;
info.reason = matrixReason(assembleMatrix, params, p.index.nTotal);
info.deviceSize = size(Vxy);
info.materialType = readField(mat, 'type', 'unknown');
info.relativeCoordinates = p.relativeCoordinates;
info.dofOrder = p.index.order;
info.totalDof = p.index.nTotal;
info.centerDof = p.dg.nCenterDof;
info.relativeDof = p.relative.nDof;
info.chunkSize = data.chunkSize;
info.potentialGrid = potentialInfo;
info.potentialDiscretization = potentialMode;
info.potentialDifference = ...
    'DeltaV = V(X+rho_x/2,Y+rho_y/2) - V(X-rho_x/2,Y-rho_y/2).';
info.CAP = cap;
info.CAPMatrixSize = size(cap.Crho);
info.CAPNnz = nnz(cap.Crho);
info.CAPNote = ['CAP remains a diagonal relative-coordinate absorber and ', ...
    'is not part of any physical X/Y boundary flux.'];
info.scale.Qdrift = data.Qdrift;
info.scale.driftScale = data.driftScale;
info.scale.capScale = data.capScale;
info.scale.potentialCoefficient = data.potentialCoefficient;
info.scale.capCoefficient = data.capCoefficient;
info.scale.note = ['Full-2D uses full2D_drift_scale only. Legacy ', ...
    'rho_drift_scale belongs to the 1D rho solver and is intentionally ', ...
    'not inherited here.'];
if hasParam(params, 'rho_drift_scale') && ~hasParam(params, 'full2D_drift_scale')
    info.scale.legacyRhoDriftScaleIgnored = params.rho_drift_scale;
end

switch potentialMode
    case 'collocated'
        info.diagonalStored = ~isempty(data.storedDiagonal);
        info.potentialDataStored = info.diagonalStored;
        info.potentialDiscretizationNote = ['Legacy pointwise multiplication ', ...
            'at relative FV cell centers. This path is retained only for ', ...
            'comparison with earlier Full-2D runs.'];
    case 'fv-consistent'
        info.diagonalStored = false;
        info.potentialDataStored = ~isempty(data.storedPotentialDifference);
        info.relativeStencil = 'tensor-product 3-by-3 FV neighborhood';
        info.relativeVertexCount = data.nRelativeVertices;
        info.potentialQuadratureSize = data.nCenterQuadrature;
        info.centerElementCount = data.nCenterElements;
        info.nPotentialSamples = data.nPotentialSamples;
        info.relativeNormalization = data.relativeNormalization;
        info.potentialDiscretizationNote = ['FV vertex projection ', ...
            '(1/16)*P''*diag(DeltaV)*P followed by tensor-product X/Y DG ', ...
            'quadrature. The assembled and matrix-free paths implement ', ...
            'the same operator.'];
end

info.apply = @(u) applyDrift(data, u);
info.makeGPUApply = @() makeGPUDriftApply(data);
info.getDiagonal = @() getDriftDiagonal(data);
info.getRowAbsSum = @() getDriftRowAbsSum(data);
info.getRelativeBlockData = @() makeDriftRelativeBlockData(data);
info.assemble = @() assembleDriftMatrix(data);
end

function data = driftData(mat, p, params, cap, Vfield, gridX, gridY, mode)
constants = physicalConstants;
Qdefault = constants.q/constants.hbar;

data = struct;
data.mode = mode;
data.nCenter = p.dg.nCenterDof;
data.nRelative = p.relative.nDof;
data.nTotal = p.index.nTotal;
data.arraySize = p.index.arraySize;
data.centerX = repmat(p.dg.X.nodes(:), p.dg.Y.nDof, 1);
data.centerY = kron(p.dg.Y.nodes(:), ones(p.dg.X.nDof, 1));
data.relativeX = repmat(p.relative.rhoX.cells(:), p.relative.NrhoY, 1);
data.relativeY = kron(p.relative.rhoY.cells(:), ones(p.relative.NrhoX, 1));
data.capProfile = full(diag(cap.Crho));
data.potential = griddedInterpolant({gridX(:), gridY(:)}, Vfield, ...
    'linear', 'nearest');
data.potentialXRange = [gridX(1), gridX(end)];
data.potentialYRange = [gridY(1), gridY(end)];

data.Qdrift = readParam(params, 'full2D_Q_drift', Qdefault);
data.driftScale = readParam(params, 'full2D_drift_scale', 1);
data.capScale = readParam(params, 'full2D_cap_scale', 1);

% With B = DeltaV - 1i*CAP, multiplication by 1i*Q gives an imaginary
% potential term and a real damping term. This is the sign convention of
% the established 1D rho-basis implementation.
data.potentialCoefficient = 1i*data.driftScale*data.Qdrift;
data.capCoefficient = data.driftScale*data.capScale*data.Qdrift;
data.chunkSize = max(1, floor(readParam(params, ...
    'full2D_driftChunkSize', min(max(1, data.nTotal), 1e6))));
data.storedDiagonal = [];
data.storedPotentialDifference = [];
data.sourceInfo = struct( ...
    'materialType', readField(mat, 'type', 'unknown'), ...
    'coordinateScale', p.coordinateScale);
end

function data = addFVConsistentData(data, p, params)
%ADDFVCONSISTENTDATA Tensor-product FV and rectangular-DG projection data.

nx = p.relative.NrhoX;
ny = p.relative.NrhoY;
Ex = spdiags(ones(nx+1, 1)*[1, 1], [-1, 0], nx+1, nx);
Ey = spdiags(ones(ny+1, 1)*[1, 1], [-1, 0], ny+1, ny);
data.relativeProjection = sparse(kron(Ey, Ex));
data.relativeNormalization = 1/16;
data.relativeVertexDegree = full(data.relativeProjection*ones(data.nRelative, 1));

[rhoVertexX, rhoVertexY] = ndgrid( ...
    p.relative.rhoX.edges, p.relative.rhoY.edges);
data.relativeVertexX = rhoVertexX(:);
data.relativeVertexY = rhoVertexY(:);
data.nRelativeVertices = numel(data.relativeVertexX);

projection = centerProjectionData(p, params);
names = fieldnames(projection);
for id = 1:numel(names)
    data.(names{id}) = projection.(names{id});
end

data.nPotentialSamples = data.nRelativeVertices ...
    * data.nCenterQuadrature * data.nCenterElements;
targetEntries = max(1, floor(readParam(params, ...
    'full2D_driftChunkTargetEntries', 2e6)));
defaultChunk = max(1, floor(targetEntries ...
    /(data.nRelativeVertices*data.nCenterQuadrature)));
data.elementChunkSize = min(data.nCenterElements, max(1, floor(readParam( ...
    params, 'full2D_driftElementChunkSize', defaultChunk))));
end

function projection = centerProjectionData(p, params)
%CENTERPROJECTIONDATA DG interpolation/projection on rectangular elements.
%
% X is the fastest local coordinate. Quadrature and local basis indices
% therefore use the tensor order (X,Y), consistent with p.dg.MXY.

nLocalX = p.dg.X.nLocal;
nLocalY = p.dg.Y.nLocal;
defaultQX = max(4, 2*nLocalX-1);
defaultQY = max(4, 2*nLocalY-1);
commonOrder = readParam(params, 'full2D_potentialQuadratureOrder', []);
if isempty(commonOrder)
    nQuadX = readParam(params, 'full2D_potentialQuadratureX', defaultQX);
    nQuadY = readParam(params, 'full2D_potentialQuadratureY', defaultQY);
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

% For U stored as [relative DOF, local center DOF], the standard strong
% multiplication operator is evaluated as
%   Uq = U*W',  Rq = C(q)*Uq,
%   Y  = Rq*diag(w)*W*M^{-T}.
centerProjection = (interpolation.*weights)*inverseMass.';

nQuadrature = numel(weights);
nLocal = nLocalX*nLocalY;
diagonalWeights = zeros(nQuadrature, nLocal);
rowAbsWeights = zeros(nQuadrature, nLocal);
for iq = 1:nQuadrature
    ell = interpolation(iq, :);
    left = inverseMass*(weights(iq)*ell.');
    diagonalWeights(iq, :) = (left.*ell.').';
    rowAbsWeights(iq, :) = (abs(left)*sum(abs(ell))).';
end

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
projection.centerDiagonalWeights = diagonalWeights;
projection.centerRowAbsWeights = rowAbsWeights;
projection.elementCenterIds = elementCenterIds;
projection.centerElementId = centerElementId;
projection.centerLocalId = centerLocalId;
projection.centerQuadratureX = centerQuadratureX;
projection.centerQuadratureY = centerQuadratureY;
projection.quadratureOrder = [nQuadX, nQuadY];
end

function y = applyDrift(data, u)
switch data.mode
    case 'collocated'
        y = applyCollocatedDrift(data, u);
    case 'fv-consistent'
        y = applyFVConsistentDrift(data, u);
end
end

function y = applyCollocatedDrift(data, u)
if numel(u) ~= data.nTotal
    error('DG:Full2D:InvalidDriftInput', ...
        'Input has %d entries, expected %d.', numel(u), data.nTotal);
end

u = u(:);
if ~isempty(data.storedDiagonal)
    y = data.storedDiagonal.*u;
    return
end

y = complex(zeros(data.nTotal, 1, 'like', u));
for first = 1:data.chunkSize:data.nTotal
    last = min(data.nTotal, first + data.chunkSize - 1);
    ids = (first:last).';
    y(ids) = collocatedDiagonalChunk(data, ids).*u(ids);
end
end

function y = applyFVConsistentDrift(data, u)
if numel(u) ~= data.nTotal
    error('DG:Full2D:InvalidDriftInput', ...
        'Input has %d entries, expected %d.', numel(u), data.nTotal);
end

U = reshape(u, data.nRelative, data.nCenter);
Y = data.capCoefficient*(data.capProfile.*U);

for first = 1:data.elementChunkSize:data.nCenterElements
    last = min(data.nCenterElements, first + data.elementChunkSize - 1);
    elements = first:last;
    centerIds = data.elementCenterIds(:, elements);
    nChunk = numel(elements);

    Ulocal = reshape(U(:, centerIds(:)), ...
        data.nRelative, data.nLocalCenter, nChunk);
    Uquadrature = pagemtimes(Ulocal, data.centerInterpolation.');
    Uquadrature = reshape(Uquadrature, data.nRelative, []);

    deltaV = potentialDifferenceForElements(data, elements);
    vertexValues = data.relativeProjection*Uquadrature;
    relativeResidual = data.relativeProjection'*(deltaV.*vertexValues);
    relativeResidual = data.potentialCoefficient ...
        * data.relativeNormalization * relativeResidual;

    relativeResidual = reshape(relativeResidual, ...
        data.nRelative, data.nCenterQuadrature, nChunk);
    Ylocal = pagemtimes(relativeResidual, data.centerProjection);
    Y(:, centerIds(:)) = Y(:, centerIds(:)) ...
        + reshape(Ylocal, data.nRelative, []);
end
y = Y(:);
end

function applyGPU = makeGPUDriftApply(data)
switch data.mode
    case 'collocated'
        gpuDiagonal = gpuArray(getDriftDiagonal(data));
        applyGPU = @(u) applyStoredDiagonal(gpuDiagonal, u, data.nTotal);
    case 'fv-consistent'
        gpuData = makeStoredGPUData(data);
        applyGPU = @(u) applyFVConsistentDrift(gpuData, u);
end
end

function gpuData = makeStoredGPUData(data)
%MAKESTOREDGPUDATA Transfer all iteration-time data to the selected GPU.
gpuData = data;
if isempty(data.storedPotentialDifference)
    potentialDifference = materializePotentialDifference(data);
else
    potentialDifference = data.storedPotentialDifference;
end
gpuData.storedPotentialDifference = gpuArray(potentialDifference);
gpuData.relativeProjection = gpuArray(data.relativeProjection);
gpuData.relativeVertexDegree = gpuArray(data.relativeVertexDegree);
gpuData.centerInterpolation = gpuArray(data.centerInterpolation);
gpuData.centerProjection = gpuArray(data.centerProjection);
gpuData.capProfile = gpuArray(data.capProfile);
% griddedInterpolant is intentionally never called in a GPU Krylov step.
gpuData.potential = [];
end

function y = applyStoredDiagonal(diagonal, u, nTotal)
if numel(u) ~= nTotal
    error('DG:Full2D:InvalidDriftInput', ...
        'Input has %d entries, expected %d.', numel(u), nTotal);
end
y = diagonal.*u(:);
end

function G = assembleDriftMatrix(data)
switch data.mode
    case 'collocated'
        diagonal = getDriftDiagonal(data);
        G = spdiags(diagonal, 0, data.nTotal, data.nTotal);
    case 'fv-consistent'
        G = assembleFVConsistentMatrix(data);
end
end

function G = assembleFVConsistentMatrix(data)
%ASSEMBLEFVCONSISTENTMATRIX Explicit sparse matrix for small diagnostics.
rows = cell(data.nCenterElements, 1);
cols = cell(data.nCenterElements, 1);
values = cell(data.nCenterElements, 1);

for element = 1:data.nCenterElements
    deltaV = potentialDifferenceForElements(data, element);
    localMatrix = sparse(data.nRelative*data.nLocalCenter, ...
        data.nRelative*data.nLocalCenter);
    for iq = 1:data.nCenterQuadrature
        Crelative = data.relativeProjection' ...
            * spdiags(deltaV(:, iq), 0, data.nRelativeVertices, ...
            data.nRelativeVertices) * data.relativeProjection;
        Crelative = data.potentialCoefficient ...
            * data.relativeNormalization * Crelative;

        ell = data.centerInterpolation(iq, :);
        % H_q maps local center trial values to strong-form test values.
        Hcenter = data.centerProjection(iq, :).'*ell;
        localMatrix = localMatrix + kron(sparse(Hcenter), Crelative);
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

capDiagonal = repmat(data.capCoefficient*data.capProfile, data.nCenter, 1);
G = sparse(vertcat(rows{:}), vertcat(cols{:}), vertcat(values{:}), ...
    data.nTotal, data.nTotal);
G = G + spdiags(capDiagonal, 0, data.nTotal, data.nTotal);
G = sparse(G);
end

function diagonal = getDriftDiagonal(data)
switch data.mode
    case 'collocated'
        if ~isempty(data.storedDiagonal)
            diagonal = data.storedDiagonal;
        else
            diagonal = materializeCollocatedDiagonal(data);
        end
    case 'fv-consistent'
        diagonal = fvConsistentDiagonal(data);
end
end

function diagonal = fvConsistentDiagonal(data)
diagonalMatrix = data.capCoefficient ...
    * repmat(data.capProfile, 1, data.nCenter);

for first = 1:data.elementChunkSize:data.nCenterElements
    last = min(data.nCenterElements, first + data.elementChunkSize - 1);
    elements = first:last;
    centerIds = data.elementCenterIds(:, elements);
    nChunk = numel(elements);
    deltaV = potentialDifferenceForElements(data, elements);

    relativeDiagonal = data.relativeNormalization ...
        * (data.relativeProjection'*deltaV);
    relativeDiagonal = reshape(relativeDiagonal, ...
        data.nRelative, data.nCenterQuadrature, nChunk);
    localDiagonal = pagemtimes(relativeDiagonal, ...
        data.centerDiagonalWeights);
    diagonalMatrix(:, centerIds(:)) = diagonalMatrix(:, centerIds(:)) ...
        + data.potentialCoefficient ...
        * reshape(localDiagonal, data.nRelative, []);
end
diagonal = diagonalMatrix(:);
end

function rowAbsSum = getDriftRowAbsSum(data)
switch data.mode
    case 'collocated'
        rowAbsSum = abs(getDriftDiagonal(data));
    case 'fv-consistent'
        rowAbsSum = fvConsistentRowAbsBound(data);
end
end

function rowAbsSum = fvConsistentRowAbsBound(data)
%FVCONSISTENTROWABSBOUND Triangle-inequality row-sum bound.
rowAbsMatrix = abs(data.capCoefficient) ...
    * repmat(abs(data.capProfile), 1, data.nCenter);

for first = 1:data.elementChunkSize:data.nCenterElements
    last = min(data.nCenterElements, first + data.elementChunkSize - 1);
    elements = first:last;
    centerIds = data.elementCenterIds(:, elements);
    nChunk = numel(elements);
    deltaV = potentialDifferenceForElements(data, elements);

    weightedVertices = abs(deltaV).*data.relativeVertexDegree;
    relativeRowBound = data.relativeNormalization ...
        * (data.relativeProjection'*weightedVertices);
    relativeRowBound = reshape(relativeRowBound, ...
        data.nRelative, data.nCenterQuadrature, nChunk);
    localRowBound = pagemtimes(relativeRowBound, ...
        data.centerRowAbsWeights);
    rowAbsMatrix(:, centerIds(:)) = rowAbsMatrix(:, centerIds(:)) ...
        + abs(data.potentialCoefficient) ...
        * reshape(localRowBound, data.nRelative, []);
end
rowAbsSum = rowAbsMatrix(:);
end

function blockData = makeDriftRelativeBlockData(data)
blockData = struct;
blockData.nCenter = data.nCenter;
blockData.nRelative = data.nRelative;
blockData.nTotal = data.nTotal;
blockData.mode = data.mode;
blockData.getBlock = @(centerId) driftRelativeBlock(data, centerId);
blockData.getCoupledBlock = @(rowId, colId) driftRelativeBlock(data, rowId, colId);
blockData.note = ['Returns the center-block diagonal of the drift/CAP ', ...
    'operator. For fv-consistent potential discretization this includes ', ...
    'the sparse 3-by-3 relative-coordinate stencil.'];
end

function block = driftRelativeBlock(data, centerId, columnId)
if nargin < 3
    columnId = centerId;
end
if centerId < 1 || centerId > data.nCenter || centerId ~= round(centerId)
    error('DG:Full2D:InvalidCenterBlock', ...
        'Center block index must be an integer in [1,%d].', data.nCenter);
end

switch data.mode
    case 'collocated'
        if centerId ~= columnId
            block = sparse(data.nRelative,data.nRelative);
            return
        end
        first = (centerId-1)*data.nRelative + 1;
        ids = first:first+data.nRelative-1;
        diagonal = collocatedDiagonalChunk(data, ids(:));
        block = spdiags(diagonal, 0, data.nRelative, data.nRelative);
    case 'fv-consistent'
        element = data.centerElementId(centerId);
        if element ~= data.centerElementId(columnId)
            block = sparse(data.nRelative,data.nRelative);
            return
        end
        localId = data.centerLocalId(centerId);
        columnLocalId = data.centerLocalId(columnId);
        deltaV = potentialDifferenceForElements(data, element);
        weights = data.centerProjection(:,localId) ...
            .*data.centerInterpolation(:,columnLocalId);
        vertexPotential = deltaV*weights;
        block = data.potentialCoefficient*data.relativeNormalization ...
            * (data.relativeProjection' ...
            * spdiags(vertexPotential, 0, data.nRelativeVertices, ...
            data.nRelativeVertices) * data.relativeProjection);
        if centerId == columnId
            block = block + spdiags(data.capCoefficient*data.capProfile, ...
                0, data.nRelative, data.nRelative);
        end
end
block = sparse(block);
end

function diagonal = materializeCollocatedDiagonal(data)
diagonal = complex(zeros(data.nTotal, 1));
for first = 1:data.chunkSize:data.nTotal
    last = min(data.nTotal, first + data.chunkSize - 1);
    ids = (first:last).';
    diagonal(ids) = collocatedDiagonalChunk(data, ids);
end
end

function diagonal = collocatedDiagonalChunk(data, ids)
relativeIds = mod(ids-1, data.nRelative) + 1;
centerIds = floor((ids-1)/data.nRelative) + 1;

rhoX = data.relativeX(relativeIds);
rhoY = data.relativeY(relativeIds);
X = data.centerX(centerIds);
Y = data.centerY(centerIds);

Vplus = evaluateExtendedPotential(data, ...
    X + 0.5*rhoX, Y + 0.5*rhoY);
Vminus = evaluateExtendedPotential(data, ...
    X - 0.5*rhoX, Y - 0.5*rhoY);
deltaV = Vplus - Vminus;

diagonal = data.potentialCoefficient*deltaV ...
    + data.capCoefficient*data.capProfile(relativeIds);
end

function deltaV = materializePotentialDifference(data)
deltaV = zeros(data.nRelativeVertices, ...
    data.nCenterQuadrature*data.nCenterElements);
for first = 1:data.elementChunkSize:data.nCenterElements
    last = min(data.nCenterElements, first + data.elementChunkSize - 1);
    elements = first:last;
    columns = elementQuadratureColumns(data, elements);
    deltaV(:, columns) = evaluatePotentialDifference(data, elements);
end
end

function deltaV = potentialDifferenceForElements(data, elements)
columns = elementQuadratureColumns(data, elements);
if ~isempty(data.storedPotentialDifference)
    deltaV = data.storedPotentialDifference(:, columns);
else
    deltaV = evaluatePotentialDifference(data, elements);
end
end

function columns = elementQuadratureColumns(data, elements)
columns = (1:data.nCenterQuadrature).' ...
    + (elements(:).'-1)*data.nCenterQuadrature;
columns = columns(:);
end

function deltaV = evaluatePotentialDifference(data, elements)
centerX = data.centerQuadratureX(:, elements);
centerY = data.centerQuadratureY(:, elements);
centerX = centerX(:).';
centerY = centerY(:).';

rhoX = 0.5*data.relativeVertexX;
rhoY = 0.5*data.relativeVertexY;
Vplus = evaluateExtendedPotential(data, ...
    centerX + rhoX, centerY + rhoY);
Vminus = evaluateExtendedPotential(data, ...
    centerX - rhoX, centerY - rhoY);
deltaV = Vplus - Vminus;
end

function values = evaluateExtendedPotential(data, queryX, queryY)
%EVALUATEEXTENDEDPOTENTIAL Constant extension in each coordinate.
%
% Clamp X and Y independently before interpolation. In particular, a
% query outside the X range retains its tangential Y coordinate instead of
% selecting the nearest two-dimensional material-grid point.
queryX = min(max(queryX, data.potentialXRange(1)), ...
    data.potentialXRange(2));
queryY = min(max(queryY, data.potentialYRange(1)), ...
    data.potentialYRange(2));
values = data.potential(queryX, queryY);
end

function [Vfield, gridX, gridY, info] = potentialOnMaterialGrid(mat, p, Vxy)
gridX = getCoordinateVector(mat, 'x', 'dx', 'Nx')*p.coordinateScale;
gridY = getCoordinateVector(mat, 'y', 'dy', 'Ny')*p.coordinateScale;
Vfield = fieldToXY(Vxy, mat);

[gridX, ix] = sort(gridX(:));
[gridY, iy] = sort(gridY(:));
Vfield = full(Vfield(ix, iy));

info = struct;
info.size = size(Vfield);
info.xRange = [gridX(1), gridX(end)];
info.yRange = [gridY(1), gridY(end)];
info.interpolation = ['linear after independent coordinate clamping ', ...
    'to the material grid'];
info.note = ['Potential values are sampled at X +/- rho_x/2 and ', ...
    'Y +/- rho_y/2. Each query coordinate is clamped independently. ', ...
    'This gives a constant continuation normal to a boundary while ', ...
    'preserving the linearly interpolated tangential boundary profile.'];
end

function Vfield = fieldToXY(field, mat)
field = squeeze(field);
if isempty(field)
    Vfield = zeros(mat.Nx, mat.Ny);
elseif isvector(field)
    if numel(field) == mat.Nx
        Vfield = repmat(field(:), 1, mat.Ny);
    elseif numel(field) == mat.Ny
        Vfield = repmat(field(:).', mat.Nx, 1);
    elseif numel(field) == mat.Nx*mat.Ny
        Vfield = reshape(field, mat.Nx, mat.Ny);
    else
        error('DG:Full2D:InvalidPotentialSize', ...
            'Vector potential has %d entries, expected Nx, Ny or Nx*Ny.', ...
            numel(field));
    end
elseif ismatrix(field)
    if isequal(size(field), [mat.Nx, mat.Ny])
        Vfield = field;
    elseif isequal(size(field), [mat.Ny, mat.Nx])
        Vfield = field.';
    elseif size(field, 1) == mat.Nx && size(field, 2) == 1
        Vfield = repmat(field, 1, mat.Ny);
    elseif size(field, 1) == 1 && size(field, 2) == mat.Ny
        Vfield = repmat(field, mat.Nx, 1);
    else
        error('DG:Full2D:InvalidPotentialSize', ...
            'Potential size is [%d, %d], expected [%d, %d].', ...
            size(field, 1), size(field, 2), mat.Nx, mat.Ny);
    end
else
    error('DG:Full2D:InvalidPotentialRank', ...
        'Full-2D drift expects a 2D potential field V(X,Y).');
end
end

function mode = normalizePotentialDiscretization(value)
mode = lower(strrep(char(value), '_', '-'));
switch mode
    case {'fv', 'fv-consistent', 'finite-volume', '1d-rho-like'}
        mode = 'fv-consistent';
    case {'collocated', 'diagonal', 'legacy'}
        mode = 'collocated';
    otherwise
        error('DG:Full2D:UnknownPotentialDiscretization', ...
            ['Unknown full2D_potentialDiscretization "%s". Use ', ...
            '"fv-consistent" or "collocated".'], char(value));
end
end

function tf = shouldAssembleDriftMatrix(params, nTotal)
if readLogicalParam(params, 'full2D_gpu', false)
    tf = false;
    return
end
mode = readParam(params, 'full2D_assembleDriftMatrix', ...
    readParam(params, 'full2D_assembleSystemMatrix', 'auto'));
tf = resolveBoolMode(mode, nTotal, readParam(params, ...
    'full2D_maxAssembledDof', 50000), 'full2D_assembleDriftMatrix');
end

function tf = shouldStoreDriftDiagonal(params, nTotal)
if readLogicalParam(params, 'full2D_gpu', false)
    tf = false;
    return
end
mode = readParam(params, 'full2D_storeDriftDiagonal', 'auto');
tf = resolveBoolMode(mode, nTotal, readParam(params, ...
    'full2D_maxStoredDriftDiagonalDof', 2e7), ...
    'full2D_storeDriftDiagonal');
end

function tf = shouldStorePotentialDifference(params, nSamples)
if readLogicalParam(params, 'full2D_gpu', false)
    % The GPU factory materializes a temporary host array and retains only
    % its device copy, avoiding a permanent duplicate in host memory.
    tf = false;
    return
end
mode = readParam(params, 'full2D_storePotentialDifference', 'auto');
tf = resolveBoolMode(mode, nSamples, readParam(params, ...
    'full2D_maxStoredPotentialDifferenceEntries', 5e7), ...
    'full2D_storePotentialDifference');
end

function tf = resolveBoolMode(mode, count, maxCount, name)
if islogical(mode)
    tf = mode;
elseif isnumeric(mode)
    tf = mode ~= 0;
else
    switch lower(char(mode))
        case 'auto'
            tf = count <= maxCount;
        case {'true', 'yes', 'on', 'assemble', 'always', 'store'}
            tf = true;
        case {'false', 'no', 'off', 'matrix-free', 'never'}
            tf = false;
        otherwise
            error('DG:Full2D:UnknownMode', ...
                'Unknown %s mode "%s".', name, char(mode));
    end
end
end

function reason = matrixReason(assembled, params, nTotal)
mode = readParam(params, 'full2D_assembleDriftMatrix', ...
    readParam(params, 'full2D_assembleSystemMatrix', 'auto'));
if assembled
    reason = sprintf('Drift matrix assembled in %s mode for %d total DOFs.', ...
        char(string(mode)), nTotal);
else
    reason = sprintf(['Drift kept matrix-free in %s mode for %d total DOFs. ', ...
        'Use info.apply(u), info.getDiagonal(), or info.assemble() for small tests.'], ...
        char(string(mode)), nTotal);
end
end

function coord = getCoordinateVector(mat, coordName, spacingName, countName)
if isfield(mat, coordName) && ~isempty(mat.(coordName))
    coord = mat.(coordName)(:);
else
    n = mat.(countName);
    d = mat.(spacingName);
    coord = (0:n-1).'*d;
end
if numel(coord) < 2
    error('DG:Full2D:InvalidGrid', '%s needs at least two grid points.', coordName);
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
if isstruct(params)
    if isfield(params, name) && ~isempty(params.(name))
        value = params.(name);
    elseif isfield(params, 'full2D') && isstruct(params.full2D) ...
            && isfield(params.full2D, name) && ~isempty(params.full2D.(name))
        value = params.full2D.(name);
    end
end
end

function value = readLogicalParam(params, name, defaultValue)
rawValue = readParam(params, name, defaultValue);
if islogical(rawValue)
    value = rawValue;
elseif isnumeric(rawValue)
    value = rawValue ~= 0;
else
    value = any(strcmpi(char(rawValue), {'true', 'yes', 'on', '1'}));
end
end

function tf = hasParam(params, name)
tf = false;
if isstruct(params)
    tf = (isfield(params, name) && ~isempty(params.(name))) ...
        || (isfield(params, 'full2D') && isstruct(params.full2D) ...
            && isfield(params.full2D, name) && ~isempty(params.full2D.(name)));
end
end

function value = readField(s, name, defaultValue)
value = defaultValue;
if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
    value = s.(name);
end
end

function c = physicalConstants
c.q = 1.602176634e-19;
c.hbar = 1.054571817e-34;
end
