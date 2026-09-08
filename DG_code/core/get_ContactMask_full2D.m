function [contactMask, info] = get_ContactMask_full2D(mat, p, sideName, yFace)
%GET_CONTACTMASK_FULL2D Mask Source/Drain reservoirs on the outer X faces.
%
% The Wigner transport domain may include oxide. Source/Drain inflow must
% therefore not be applied to the entire vertical outer rectangle face, but
% only to the semiconductor contact opening. The default mask is derived
% from mat.Nd on the X-left/X-right edge. This map is created by buildDevice
% from its geometry table, so variable widths and oxide thicknesses are
% respected without duplicating those dimensions here.
%
% Explicit masks mat.dg.params.full2D_sourceContactMask and
% mat.dg.params.full2D_drainContactMask override the material-derived mask.

params = getDGParams(mat);
explicitName = ['full2D_', lower(sideName), 'ContactMask'];
yFace = yFace(:);

if isfield(params, explicitName) && ~isempty(params.(explicitName))
    [contactMask, origin] = explicitContactMask(params.(explicitName), mat, p, yFace);
elseif hasUsableField(mat, 'Nd')
    [contactMask, origin] = maskFromDoping(mat, p, sideName, yFace, params);
else
    [contactMask, origin] = maskFromMaterialOrGeometry(mat, p, sideName, yFace);
end

if ~any(contactMask)
    [contactMask, fallbackOrigin] = maskFromMaterialOrGeometry(mat, p, sideName, yFace);
    origin = [origin, ' -> ', fallbackOrigin];
end

info = struct;
info.side = sideName;
info.origin = origin;
info.nFaceDof = numel(contactMask);
info.nContactDof = nnz(contactMask);
info.nClosedDof = nnz(~contactMask);
info.closedBoundaryType = 'specular-reflection in rho_x on masked X-face segments';
info.note = ['Mask entries true receive Source/Drain inflow. Masked-off ', ...
    'outer X-face entries are closed and use reflected ghost data instead.'];
end

function [mask, origin] = explicitContactMask(maskData, mat, p, yFace)
maskData = maskData(:);
if numel(maskData) == numel(yFace)
    mask = logical(maskData);
    origin = 'explicit-face-mask';
elseif isfield(mat, 'Ny') && numel(maskData) == mat.Ny
    gridY = getCoordinateVector(mat, 'y', 'dy', 'Ny')*readField(p, 'coordinateScale', 1);
    mask = interpLogical(maskData, gridY, yFace);
    origin = 'explicit-grid-y-mask';
else
    error('DG:Full2D:InvalidContactMask', ...
        'Contact mask has %d entries, expected %d face DOFs or mat.Ny entries.', ...
        numel(maskData), numel(yFace));
end
end

function [mask, origin] = maskFromDoping(mat, p, sideName, yFace, params)
threshold = readParam(params, 'full2D_contactNdThreshold', 0);
Nd = fieldToXY(mat.Nd, mat, 0);
sideIndex = sideToXIndex(mat, sideName);
gridY = getCoordinateVector(mat, 'y', 'dy', 'Ny')*readField(p, 'coordinateScale', 1);
mask = interpLogical(Nd(sideIndex, :) > threshold, gridY, yFace);
origin = 'buildDevice-Nd-edge';
end

function [mask, origin] = maskFromMaterialOrGeometry(mat, p, sideName, yFace)
if hasUsableField(mat, 'Eg') && hasUsableField(mat, 'Eg_ox')
    Eg = fieldToXY(mat.Eg, mat, readField(mat, 'Eg_c', 0));
    sideIndex = sideToXIndex(mat, sideName);
    gridY = getCoordinateVector(mat, 'y', 'dy', 'Ny')*readField(p, 'coordinateScale', 1);
    threshold = 0.5*(readField(mat, 'Eg_c', min(Eg(:))) + mat.Eg_ox);
    mask = interpLogical(Eg(sideIndex, :) < threshold, gridY, yFace);
    origin = 'buildDevice-material-edge';
elseif isfield(mat, 'W_ox') && isfield(mat, 'W_c')
    scale = readField(p, 'coordinateScale', 1);
    mask = yFace >= 0 & yFace <= mat.W_c*scale;
    origin = 'buildDevice-W_ox-W_c-geometry';
elseif isfield(mat, 'W_oxy') && isfield(mat, 'W_cy')
    scale = readField(p, 'coordinateScale', 1);
    mask = yFace >= 0 & yFace <= mat.W_cy*scale;
    origin = 'buildDevice-W_oxy-W_cy-geometry';
else
    mask = true(size(yFace));
    origin = 'all-face-dofs-no-material-mask-available';
end
mask = logical(mask(:));
end

function mask = interpLogical(maskGrid, gridY, yFace)
maskGrid = double(maskGrid(:));
mask = interp1(gridY(:), maskGrid, yFace, 'nearest', 'extrap') > 0.5;
mask = logical(mask(:));
end

function field = fieldToXY(rawField, mat, defaultValue)
rawField = squeeze(rawField);
fieldSize = size(rawField);

if isempty(rawField)
    field = defaultValue*ones(mat.Nx, mat.Ny);
elseif ismatrix(rawField) && isequal(fieldSize, [mat.Nx, mat.Ny])
    field = rawField;
elseif ismatrix(rawField) && isequal(fieldSize, [mat.Ny, mat.Nx])
    field = rawField.';
elseif isvector(rawField) && numel(rawField) == mat.Nx
    field = repmat(rawField(:), 1, mat.Ny);
elseif isvector(rawField) && numel(rawField) == mat.Ny
    field = repmat(rawField(:).', mat.Nx, 1);
elseif ndims(rawField) == 3 && size(rawField, 2) == mat.Nx && size(rawField, 3) == mat.Ny
    field = squeeze(rawField(1, :, :));
elseif ndims(rawField) == 4 && size(rawField, 2) == mat.Nx && size(rawField, 3) == mat.Ny
    zIndex = ceil(size(rawField, 4)/2);
    field = squeeze(rawField(1, :, :, zIndex));
else
    error('DG:Full2D:InvalidMaterialField', ...
        'Cannot interpret material field with size [%s] as an X/Y field.', ...
        num2str(size(rawField)));
end
end

function sideIndex = sideToXIndex(mat, sideName)
switch lower(sideName)
    case 'source'
        sideIndex = 1;
    case 'drain'
        sideIndex = mat.Nx;
    otherwise
        error('DG:Full2D:UnknownReservoirSide', ...
            'Unknown reservoir side "%s". Use "source" or "drain".', sideName);
end
end

function coord = getCoordinateVector(mat, coordName, spacingName, countName)
if isfield(mat, coordName) && ~isempty(mat.(coordName))
    coord = mat.(coordName)(:).';
else
    coord = (0:mat.(countName)-1)*mat.(spacingName);
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

function tf = hasUsableField(s, name)
tf = isstruct(s) && isfield(s, name) && ~isempty(s.(name));
end
