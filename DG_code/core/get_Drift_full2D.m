function [G, info] = get_Drift_full2D(mat, p, Vxy)
%GET_DRIFT_FULL2D Local potential and CAP operator for Full-2D Wigner-DG.
%
% Unknown and vector order:
%   F(iX, iY, iRhoX, iRhoY) is stored as F(:). X-DG is the fastest index,
%   followed by Y-DG, rho_x FV cells and rho_y FV cells.
%
% In the rho-basis prototype used here, the nonlocal Wigner potential is
% represented as a collocated relative-coordinate potential difference
%
%   DeltaV = V(X+rho_x/2,Y+rho_y/2) - V(X-rho_x/2,Y-rho_y/2).
%
% This gives a diagonal global operator in the current discrete basis. The
% CAP is added separately as a relative-coordinate absorber:
%
%   C_rho = C_x kron I_y + I_x kron C_y.
%
% No physical Source/Drain or specular-wall boundary contribution is mixed
% into this routine. Those flux terms live in get_Diff_full2D.

if nargin < 3 || isempty(Vxy)
    if isfield(mat, 'V') && ~isempty(mat.V)
        Vxy = mat.V;
    else
        Vxy = zeros(mat.Nx, mat.Ny);
    end
end

params = getDGParams(mat);
cap = get_CAP_full2D(p);
[Vfield, gridX, gridY, potentialInfo] = potentialOnMaterialGrid(mat, p, Vxy);

data = driftData(mat, p, params, cap, Vfield, gridX, gridY);
assembleMatrix = shouldAssembleDriftMatrix(params, p.index.nTotal);
storeDiagonal = assembleMatrix || shouldStoreDriftDiagonal(params, p.index.nTotal);

if storeDiagonal
    data.storedDiagonal = materializeDriftDiagonal(data);
else
    data.storedDiagonal = [];
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
info.diagonalStored = ~isempty(data.storedDiagonal);
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
info.potentialDifference = ...
    'DeltaV = V(X+rho_x/2,Y+rho_y/2) - V(X-rho_x/2,Y-rho_y/2).';
info.CAP = cap;
info.CAPMatrixSize = size(cap.Crho);
info.CAPNnz = nnz(cap.Crho);
info.CAPNote = ['CAP is separable in rho_x/rho_y and remains separate ', ...
    'from physical Source/Drain and specular-reflection fluxes.'];
info.scale.Qdrift = data.Qdrift;
info.scale.driftScale = data.driftScale;
info.scale.capScale = data.capScale;
info.scale.potentialCoefficient = data.potentialCoefficient;
info.scale.capCoefficient = data.capCoefficient;
info.apply = @(u) applyDriftDiagonal(data, u);
info.getDiagonal = @() getDriftDiagonal(data);
info.assemble = @() assembleDriftMatrix(data);
end

function data = driftData(mat, p, params, cap, Vfield, gridX, gridY)
constants = physicalConstants;
Qdefault = constants.q/constants.hbar;

data = struct;
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

data.Qdrift = readParam(params, 'full2D_Q_drift', Qdefault);
data.driftScale = readParam(params, 'full2D_drift_scale', ...
    readParam(params, 'rho_drift_scale', 1));
data.capScale = readParam(params, 'full2D_cap_scale', 1);

% With B = DeltaV - 1i*CAP, the Wigner drift contribution is 1i*Q*B.
% Therefore the potential part is imaginary and the absorber enters as a
% real damping contribution. This mirrors the sign convention used by the
% existing 1D rho-basis drift implementation.
data.potentialCoefficient = 1i*data.driftScale*data.Qdrift;
data.capCoefficient = data.driftScale*data.capScale*data.Qdrift;
data.chunkSize = max(1, floor(readParam(params, ...
    'full2D_driftChunkSize', min(max(1, data.nTotal), 1e6))));
data.storedDiagonal = [];
data.sourceInfo = struct( ...
    'materialType', readField(mat, 'type', 'unknown'), ...
    'coordinateScale', p.coordinateScale);
end

function y = applyDriftDiagonal(data, u)
if numel(u) ~= data.nTotal
    error('DG:Full2D:InvalidDriftInput', ...
        'Input has %d entries, expected %d.', numel(u), data.nTotal);
end

u = u(:);
if ~isempty(data.storedDiagonal)
    y = data.storedDiagonal.*u;
    return
end

y = complex(zeros(data.nTotal, 1));
for first = 1:data.chunkSize:data.nTotal
    last = min(data.nTotal, first + data.chunkSize - 1);
    ids = (first:last).';
    y(ids) = driftDiagonalChunk(data, ids).*u(ids);
end
end

function G = assembleDriftMatrix(data)
diagonal = data.storedDiagonal;
if isempty(diagonal)
    diagonal = materializeDriftDiagonal(data);
end
G = spdiags(diagonal, 0, data.nTotal, data.nTotal);
end

function diagonal = getDriftDiagonal(data)
if ~isempty(data.storedDiagonal)
    diagonal = data.storedDiagonal;
else
    diagonal = materializeDriftDiagonal(data);
end
end

function diagonal = materializeDriftDiagonal(data)
diagonal = complex(zeros(data.nTotal, 1));
for first = 1:data.chunkSize:data.nTotal
    last = min(data.nTotal, first + data.chunkSize - 1);
    ids = (first:last).';
    diagonal(ids) = driftDiagonalChunk(data, ids);
end
end

function diagonal = driftDiagonalChunk(data, ids)
%DRIFTDIAGONALCHUNK Evaluate the local diagonal entries for linear DOFs.
%
% The linear index is split into a center index and a relative index using
% the global ordering (center fastest, relative slowest).

centerIds = mod(ids-1, data.nCenter) + 1;
relativeIds = floor((ids-1)/data.nCenter) + 1;

rhoX = data.relativeX(relativeIds);
rhoY = data.relativeY(relativeIds);
X = data.centerX(centerIds);
Y = data.centerY(centerIds);

Vplus = data.potential(X + 0.5*rhoX, Y + 0.5*rhoY);
Vminus = data.potential(X - 0.5*rhoX, Y - 0.5*rhoY);
deltaV = Vplus - Vminus;

diagonal = data.potentialCoefficient*deltaV ...
    + data.capCoefficient*data.capProfile(relativeIds);
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
info.interpolation = 'linear inside the material grid, nearest outside';
info.note = ['Potential values are sampled at X +/- rho_x/2 and ', ...
    'Y +/- rho_y/2 for the local rho-basis prototype.'];
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
        'Full-2D drift currently expects a 2D potential field V(X,Y).');
end
end

function tf = shouldAssembleDriftMatrix(params, nTotal)
mode = readParam(params, 'full2D_assembleDriftMatrix', ...
    readParam(params, 'full2D_assembleSystemMatrix', 'auto'));
tf = resolveBoolMode(mode, nTotal, readParam(params, ...
    'full2D_maxAssembledDof', 50000), 'full2D_assembleDriftMatrix');
end

function tf = shouldStoreDriftDiagonal(params, nTotal)
mode = readParam(params, 'full2D_storeDriftDiagonal', 'auto');
tf = resolveBoolMode(mode, nTotal, readParam(params, ...
    'full2D_maxStoredDriftDiagonalDof', 2e7), ...
    'full2D_storeDriftDiagonal');
end

function tf = resolveBoolMode(mode, nTotal, maxDof, name)
if islogical(mode)
    tf = mode;
elseif isnumeric(mode)
    tf = mode ~= 0;
else
    switch lower(char(mode))
        case 'auto'
            tf = nTotal <= maxDof;
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
