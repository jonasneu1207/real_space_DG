function [rhoBoundary, info] = get_DefaultReservoirRho_full2D(mat, p, sideName, Ef, Vxy)
%GET_DEFAULTRESERVOIRRHO_FULL2D Material-based Source/Drain data in rho.
%
% The full-2D unknown is ordered as (rho_x, rho_y), with rho_x fastest. This
% helper builds a contact Fermi distribution on a separable (k_x,k_y) grid
% and transforms it to the relative-coordinate basis by the discrete
% Fourier-cosine quadrature
%
%   rho(rho_x,rho_y) = int int f(k_x,k_y)
%       cos(rho_x*k_x) cos(rho_y*k_y) dk_x dk_y/(2*pi)^2.
%
% Source and Drain are evaluated at the physical X-left/X-right contact. The
% contact profiles V, me_x and me_y are sampled along all Y-DG face DOFs, so
% the returned matrix has one rho-boundary column per Y boundary node.

if nargin < 5 || isempty(Vxy)
    Vxy = readField(mat, 'V', []);
end

params = getDGParams(mat);
if ~readParam(params, 'full2D_useMaterialReservoir', true)
    [rhoBoundary, info] = zeroReservoir(p, sideName, ...
        'full2D_useMaterialReservoir is false.');
    return
end

if isempty(Ef) || any(~isfinite(Ef(:)))
    [rhoBoundary, info] = zeroReservoir(p, sideName, ...
        'No finite reservoir Fermi level was supplied.');
    return
end

constants = physicalConstants;
temp = readField(mat, 'Temp', 300);
degFactor = readDegeneracy(mat);

nKx = readParam(params, 'full2D_reservoir_Nkx', p.relative.NrhoX);
nKy = readParam(params, 'full2D_reservoir_Nky', p.relative.NrhoY);
[kx, dkx] = makeKGrid(nKx, p.relative.rhoX.length);
[ky, dky] = makeKGrid(nKy, p.relative.rhoY.length);

rhoX = p.relative.rhoX.cells(:);
rhoY = p.relative.rhoY.cells(:);
cosX = cos(rhoX*kx);
cosY = cos(rhoY*ky);

yFace = p.dg.Y.nodes(:);
[Vedge, mXedge, mYedge, profileInfo] = contactProfiles(mat, p, Vxy, sideName, yFace);

rhoBoundary = zeros(p.relative.nDof, numel(yFace));
for iy = 1:numel(yFace)
    rhoGrid = transformReservoirColumn(cosX, cosY, kx, ky, dkx, dky, ...
        Vedge(iy), mXedge(iy), mYedge(iy), Ef, temp, degFactor, constants);
    rhoBoundary(:, iy) = rhoGrid(:);
end

info = struct;
info.origin = 'material-default';
info.side = sideName;
info.basis = 'rho_x/rho_y';
info.relativeOrder = {'rho_x', 'rho_y'};
info.note = ['Default full-2D Source/Drain data were generated from ', ...
    'Ef, contact potential, effective masses, temperature and degeneracy.'];
info.transform = 'separable cosine Fourier quadrature from k_x/k_y to rho_x/rho_y';
info.kxRange = [min(kx), max(kx)];
info.kyRange = [min(ky), max(ky)];
info.dk = [dkx, dky];
info.temperature = temp;
info.degeneracy = degFactor;
info.contactProfiles = profileInfo;
info.size = size(rhoBoundary);
end

function rhoGrid = transformReservoirColumn(cosX, cosY, kx, ky, dkx, dky, ...
        Vedge, mX, mY, Ef, temp, degFactor, constants)
[KX, KY] = ndgrid(kx, ky);
EkX = constants.hbar^2*KX.^2/(2*mX*constants.q);
EkY = constants.hbar^2*KY.^2/(2*mY*constants.q);
energy = Vedge + EkX + EkY;

arg = constants.q*(energy - Ef)/(constants.kB*temp);
arg = min(max(arg, -80), 80);
phaseDistribution = 2*degFactor./(1 + exp(arg));

rhoGrid = cosX*phaseDistribution*cosY.'*dkx*dky/(2*pi)^2;
end

function [Vedge, mXedge, mYedge, info] = contactProfiles(mat, p, Vxy, sideName, yFace)
sideIndex = sideToXIndex(mat, sideName);
gridY = getCoordinateVector(mat, 'y', 'dy', 'Ny')*readField(p, 'coordinateScale', 1);

Vfield = fieldToXY(Vxy, mat, 0);
mXfield = fieldToXY(readField(mat, 'me_x', []), mat, readField(mat, 'me_x_ch', 0.041));
mYfield = fieldToXY(readField(mat, 'me_y', []), mat, readField(mat, 'me_y_ch', 0.041));

Vedge = interpolateEdge(Vfield(sideIndex, :), gridY, yFace);
mXedge = massToKg(interpolateEdge(mXfield(sideIndex, :), gridY, yFace));
mYedge = massToKg(interpolateEdge(mYfield(sideIndex, :), gridY, yFace));

info = struct;
info.side = sideName;
info.xIndex = sideIndex;
info.yFaceCount = numel(yFace);
info.VRange = [min(Vedge), max(Vedge)];
info.meXRangeKg = [min(mXedge), max(mXedge)];
info.meYRangeKg = [min(mYedge), max(mYedge)];
info.axis = 'Y-DG face nodes along the Source/Drain contact';
end

function data = interpolateEdge(edgeData, gridY, yFace)
edgeData = edgeData(:);
if isscalar(edgeData)
    data = repmat(edgeData, size(yFace));
else
    data = interp1(gridY(:), edgeData, yFace, 'linear', 'extrap');
end
data = data(:);
end

function field = fieldToXY(rawField, mat, defaultValue)
if isempty(rawField)
    field = defaultValue*ones(mat.Nx, mat.Ny);
    return
end

rawField = squeeze(rawField);
fieldSize = size(rawField);

if ismatrix(rawField) && isequal(fieldSize, [mat.Nx, mat.Ny])
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

function massKg = massToKg(massValue)
constants = physicalConstants;
massKg = massValue;
relativeMask = abs(massValue) > 1e-25;
massKg(relativeMask) = massValue(relativeMask)*constants.m0;
massKg(massKg <= 0) = 0.041*constants.m0;
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

function [k, dk] = makeKGrid(nK, lengthValue)
if nK < 2
    error('DG:Full2D:InvalidKGrid', 'Reservoir k grid needs at least two points.');
end
k = -2*pi/lengthValue*((1:nK) - 0.5*(nK+1));
dk = abs(k(2)-k(1));
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

function degFactor = readDegeneracy(mat)
degFactor = 1;
if isfield(mat, 'deg_factor') && ~isempty(mat.deg_factor)
    degFactor = mat.deg_factor(1);
end
end

function c = physicalConstants
c.hbar = 1.054571817e-34;
c.q = 1.602176634e-19;
c.kB = 1.38064852e-23;
c.m0 = 9.1093837015e-31;
end

function [rhoBoundary, info] = zeroReservoir(p, sideName, note)
rhoBoundary = zeros(p.relative.nDof, p.dg.Y.nDof);
info = struct;
info.origin = 'zero-placeholder';
info.side = sideName;
info.note = note;
info.size = size(rhoBoundary);
end
