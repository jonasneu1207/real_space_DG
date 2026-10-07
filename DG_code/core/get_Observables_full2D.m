function [n, jx, jy, info] = get_Observables_full2D(rho, p, mat)
%GET_OBSERVABLES_FULL2D Density and current diagnostics for Full-2D rho data.
%
% The full unknown is rho(rho_x,rho_y,X,Y), vectorized with rho_x fastest.
% The carrier density is sampled at zero relative separation:
%
%   n(X,Y) = real(rho(0,0,X,Y)).
%
% The current densities are estimated from the derivative of the density
% matrix at zero relative separation, analogous to the existing 1D rho
% solver:
%
%   jx proportional imag(d rho / d rho_x)|_0,
%   jy proportional imag(d rho / d rho_y)|_0.
%
% The returned n, jx and jy are projected onto the material Y grid when that
% grid is available. Raw rectangular DG values are kept in info.*DG.

if numel(rho) ~= p.index.nTotal
    error('DG:Full2D:InvalidObservableInput', ...
        'rho has %d entries, expected %d.', numel(rho), p.index.nTotal);
end

rho = rho(:);
rhoByCenter = reshape(rho, p.relative.nDof, p.dg.nCenterDof);
rho4D = reshape(rho, p.index.arraySize);

[rhoXIds, rhoXWeights, rhoXMode] = zeroInterpolationWeights(p.relative.rhoX.cells);
[rhoYIds, rhoYWeights, rhoYMode] = zeroInterpolationWeights(p.relative.rhoY.cells);
rhoXDerivativeRow = derivativeRowAtZero(p.relative.rhoX.cells, p.relative.DrhoX);
rhoYDerivativeRow = derivativeRowAtZero(p.relative.rhoY.cells, p.relative.DrhoY);

nDG = zeros(p.dg.X.nDof, p.dg.Y.nDof);
for ix = 1:numel(rhoXIds)
    for iy = 1:numel(rhoYIds)
        nDG = nDG + rhoXWeights(ix)*rhoYWeights(iy) ...
            *real(squeeze(rho4D(rhoXIds(ix), rhoYIds(iy), :, :)));
    end
end

derivativeX = zeros(p.dg.nCenterDof, 1);
for iy = 1:numel(rhoYIds)
    relIds = (rhoYIds(iy)-1)*p.relative.NrhoX + (1:p.relative.NrhoX);
    derivativeX = derivativeX ...
        + rhoYWeights(iy)*(rhoXDerivativeRow*rhoByCenter(relIds, :)).';
end

derivativeY = zeros(p.dg.nCenterDof, 1);
for ix = 1:numel(rhoXIds)
    relIds = rhoXIds(ix) + (0:p.relative.NrhoY-1)*p.relative.NrhoX;
    derivativeY = derivativeY ...
        + rhoXWeights(ix)*(rhoYDerivativeRow*rhoByCenter(relIds, :)).';
end

[scaleX, scaleY, currentScaleInfo] = currentScales(mat, p);
jxDG = scaleX.*imag(reshape(derivativeX, p.dg.X.nDof, p.dg.Y.nDof));
jyDG = scaleY.*imag(reshape(derivativeY, p.dg.X.nDof, p.dg.Y.nDof));

[n, targetY, projectionInfo] = projectToMaterialY(nDG, p, mat);
jx = projectValuesToY(jxDG, p.dg.Y.nodes(:), targetY);
jy = projectValuesToY(jyDG, p.dg.Y.nodes(:), targetY);

info = struct;
info.currentComponents = {'jx', 'jy'};
info.nDG = nDG;
info.jxDG = jxDG;
info.jyDG = jyDG;
info.yDG = p.dg.Y.nodes(:);
info.y = targetY(:);
info.relativeZero.rhoXIds = rhoXIds;
info.relativeZero.rhoYIds = rhoYIds;
info.relativeZero.rhoXWeights = rhoXWeights;
info.relativeZero.rhoYWeights = rhoYWeights;
info.relativeZero.rhoXMode = rhoXMode;
info.relativeZero.rhoYMode = rhoYMode;
info.currentScale = currentScaleInfo;
info.projection = projectionInfo;
info.jxIntegratedOverY = integrateOverY(jx, targetY);
end

function [ids, weights, mode] = zeroInterpolationWeights(cells)
cells = cells(:);
[distance, nearest] = min(abs(cells));
tol = max(1e-14, 1e-12*max(1, max(abs(cells))));
if distance <= tol
    ids = nearest;
    weights = 1;
    mode = 'exact-cell';
    return
end

left = find(cells < 0, 1, 'last');
right = find(cells > 0, 1, 'first');
if ~isempty(left) && ~isempty(right)
    h = cells(right) - cells(left);
    ids = [left, right];
    weights = [cells(right)/h, -cells(left)/h];
    mode = 'linear-bracket';
else
    ids = nearest;
    weights = 1;
    mode = 'nearest-cell';
end
end

function row = derivativeRowAtZero(cells, D)
cells = cells(:);
[distance, nearest] = min(abs(cells));
tol = max(1e-14, 1e-12*max(1, max(abs(cells))));
if distance <= tol
    row = full(D(nearest, :));
    return
end

left = find(cells < 0, 1, 'last');
right = find(cells > 0, 1, 'first');
row = zeros(1, numel(cells));
if ~isempty(left) && ~isempty(right)
    h = cells(right) - cells(left);
    row(left) = -1/h;
    row(right) = 1/h;
elseif nearest == 1
    h = cells(2) - cells(1);
    row(1) = -1/h;
    row(2) = 1/h;
elseif nearest == numel(cells)
    h = cells(end) - cells(end-1);
    row(end-1) = -1/h;
    row(end) = 1/h;
else
    h = cells(nearest+1) - cells(nearest-1);
    row(nearest-1) = -1/h;
    row(nearest+1) = 1/h;
end
end

function [valuesOnTarget, targetY, info] = projectToMaterialY(valuesDG, p, mat)
if isfield(mat, 'y') && ~isempty(mat.y)
    targetY = mat.y(:)*p.coordinateScale;
elseif isfield(mat, 'Ny') && isfield(mat, 'dy')
    targetY = (0:mat.Ny-1).'*mat.dy*p.coordinateScale;
else
    targetY = p.dg.Y.nodes(:);
end

valuesOnTarget = projectValuesToY(valuesDG, p.dg.Y.nodes(:), targetY);

info = struct;
info.method = 'duplicate DG Y nodes are averaged, then linearly interpolated';
info.target = 'material Y grid when available';
info.nDGColumns = size(valuesDG, 2);
info.nTargetColumns = numel(targetY);
end

function projected = projectValuesToY(valuesDG, yDG, targetY)
yDG = yDG(:);
targetY = targetY(:);
[uniqueY, ~, group] = unique(yDG);
averaged = zeros(size(valuesDG, 1), numel(uniqueY));
for iy = 1:numel(uniqueY)
    averaged(:, iy) = mean(valuesDG(:, group == iy), 2);
end

if isscalar(uniqueY)
    projected = repmat(averaged, 1, numel(targetY));
else
    projected = interp1(uniqueY, averaged.', targetY, 'linear', 'extrap').';
end
end

function jIntegrated = integrateOverY(jx, y)
if isempty(y) || isscalar(y)
    jIntegrated = jx(:, 1);
else
    jIntegrated = trapz(y(:), jx, 2);
end
end

function [scaleX, scaleY, info] = currentScales(mat, p)
constants = physicalConstants;
params = getDGParams(mat);
massModel = lower(strrep(char(readParam(params, ...
    'full2D_massModel', 'constant')), '_', '-'));

if any(strcmp(massModel, {'position-dependent', 'spatial', 'variable', ...
        'variable-mass', 'bdd', 'ben-daniel-duke', ...
        'position-dependent-bdd'}))
    inverseMassX = inverseMassOnDGGrid(mat, p, 'x', ...
        readField(mat, 'me_x_ch', 0.041));
    inverseMassY = inverseMassOnDGGrid(mat, p, 'y', ...
        readField(mat, 'me_y_ch', readField(mat, 'me_x_ch', 0.041)));
    scaleX = constants.q*constants.hbar/constants.m0*inverseMassX;
    scaleY = constants.q*constants.hbar/constants.m0*inverseMassY;
    resolvedModel = 'position-dependent-bdd';
else
    massX = massToKg(readField(mat, 'me_x_ch', 0.041));
    massY = massToKg(readField(mat, 'me_y_ch', ...
        readField(mat, 'me_x_ch', 0.041)));
    scaleX = constants.q*constants.hbar/massX;
    scaleY = constants.q*constants.hbar/massY;
    inverseMassX = constants.m0/massX;
    inverseMassY = constants.m0/massY;
    resolvedModel = 'constant';
end

info = struct;
info.massModel = resolvedModel;
info.scaleXRange = [min(scaleX, [], 'all'), max(scaleX, [], 'all')];
info.scaleYRange = [min(scaleY, [], 'all'), max(scaleY, [], 'all')];
info.inverseRelativeMassXRange = [min(inverseMassX, [], 'all'), ...
    max(inverseMassX, [], 'all')];
info.inverseRelativeMassYRange = [min(inverseMassY, [], 'all'), ...
    max(inverseMassY, [], 'all')];
info.note = ['For position-dependent masses, j_x and j_y use ', ...
    'q*hbar/m_x(X,Y) and q*hbar/m_y(X,Y) on the DG center grid.'];
end

function inverseMass = inverseMassOnDGGrid(mat, p, axisName, defaultMass)
fieldName = ['me_', axisName];
mass = massFieldToXY(readField(mat, fieldName, defaultMass), ...
    mat, defaultMass);
constants = physicalConstants;
kgMask = abs(mass) < 1e-25;
mass(kgMask) = mass(kgMask)/constants.m0;
if any(~isfinite(mass), 'all') || any(mass <= 0, 'all')
    error('DG:Full2D:InvalidObservableMass', ...
        '%s must contain finite positive masses.', fieldName);
end

gridX = coordinateVector(mat, 'x', 'dx', 'Nx')*p.coordinateScale;
gridY = coordinateVector(mat, 'y', 'dy', 'Ny')*p.coordinateScale;
[gridX, ix] = sort(gridX(:));
[gridY, iy] = sort(gridY(:));
inverseField = 1./mass(ix, iy);
interpolant = griddedInterpolant({gridX, gridY}, inverseField, ...
    'linear', 'nearest');
[X, Y] = ndgrid(p.dg.X.nodes(:), p.dg.Y.nodes(:));
inverseMass = interpolant(X, Y);
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
        error('DG:Full2D:InvalidObservableMassSize', ...
            'Mass array must end in dimensions Nx-by-Ny.');
    end
    field = reshape(field, [], mat.Nx, mat.Ny);
    field = squeeze(field(1, :, :));
else
    field = squeeze(field);
end
if isequal(size(field), [mat.Nx, mat.Ny])
    % Already in the material-grid convention field(X,Y).
elseif isequal(size(field), [mat.Ny, mat.Nx])
    field = field.';
else
    error('DG:Full2D:InvalidObservableMassSize', ...
        'Mass field must have size Nx-by-Ny.');
end
end

function coord = coordinateVector(mat, name, spacingName, countName)
if isfield(mat, name) && ~isempty(mat.(name))
    coord = mat.(name)(:);
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

function value = readField(s, name, defaultValue)
value = defaultValue;
if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
    value = s.(name);
end
end

function c = physicalConstants
c.q = 1.602176634e-19;
c.hbar = 1.054571817e-34;
c.m0 = 9.1093837015e-31;
end
