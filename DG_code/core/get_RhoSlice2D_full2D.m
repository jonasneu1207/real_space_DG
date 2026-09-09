function [slice, axes, info] = get_RhoSlice2D_full2D(rho, p, freeDims, fixed)
%GET_RHOSLICE2D_FULL2D Extract a 2D slice from full-2D rho data.
%
% The Full-2D unknown is
%
%   rho = rho(rho_x, rho_y, X, Y)
%
% and the global vector order is fixed throughout the Full-2D path as
%
%   rho(iRhoX, iRhoY, iX, iY) -> rho(:),
%
% with rho_x FV as the fastest index, followed by rho_y FV, X-DG and Y-DG. This
% helper reshapes the vector with p.index.arraySize, keeps exactly two axes
% free, and fixes the other two axes by index or nearest coordinate. The
% extraction indexes the 2D slice before any permutation, so it does not
% duplicate the complete 4D solution array.
%
% Examples:
%   [S, ax] = get_RhoSlice2D_full2D(DG.rho, DG.p, ...
%       {'X', 'rho_x'}, struct('Y', 0, 'rho_y', 0));
%   imagesc(ax.rho_x, ax.X, real(S));
%
%   [S, ax] = DG.getRhoSlice2D({'rho_x', 'rho_y'}, ...
%       struct('X', 20e-9, 'Y', 0));
%
% The returned matrix dimension follows the requested freeDims order:
% size(slice) = [numel(axes.(freeDims{1})), numel(axes.(freeDims{2}))].

if nargin < 4 || isempty(fixed)
    fixed = struct;
end

if isempty(rho)
    error('DG:Full2D:EmptyRhoSolution', ...
        'rho is empty. Run a Full-2D solve before extracting a slice.');
end
if numel(rho) ~= p.index.nTotal
    error('DG:Full2D:InvalidSliceInput', ...
        'rho has %d entries, expected %d.', numel(rho), p.index.nTotal);
end

dim = dimensionInfo(p);
freeIds = parseFreeDimensions(freeDims, dim);
fixedIds = setdiff(1:4, freeIds, 'stable');

fixedIndex = ones(1, 4);
fixedValue = nan(1, 4);
fixedSource = cell(1, 4);
for id = fixedIds
    [fixedIndex(id), fixedValue(id), fixedSource{id}] = ...
        resolveFixedIndex(fixed, dim(id));
end

% Select first, then permute only the small 2D slice. That keeps the output
% matrix in exactly the order requested by freeDims, even for combinations
% such as {'rho_y','X'}, without duplicating the complete 4D array.
rho4D = reshape(rho(:), p.index.arraySize);
subs = repmat({':'}, 1, 4);
for id = fixedIds
    subs{id} = fixedIndex(id);
end

rawSlice = rho4D(subs{:});
rawShape = ones(1, 4);
rawShape(freeIds(1)) = numel(dim(freeIds(1)).coord);
rawShape(freeIds(2)) = numel(dim(freeIds(2)).coord);
rawSlice = reshape(rawSlice, rawShape);

slice = permute(rawSlice, [freeIds(:).', fixedIds(:).']);
slice = reshape(slice, numel(dim(freeIds(1)).coord), ...
    numel(dim(freeIds(2)).coord));

axes = struct;
for k = 1:2
    axes.(dim(freeIds(k)).name) = dim(freeIds(k)).coord(:);
end
axes.freeDims = {dim(freeIds(1)).name, dim(freeIds(2)).name};

info = struct;
info.freeDims = axes.freeDims;
info.fixedDims = {dim(fixedIds(1)).name, dim(fixedIds(2)).name};
info.dofOrder = p.index.order;
info.vectorOrder = ...
    'rho(iRhoX,iRhoY,iX,iY) is vectorized with rho_x fastest.';
info.indices = struct;
for id = 1:4
    field = dim(id).name;
    if ismember(id, freeIds)
        info.indices.(field) = struct( ...
            'index', [], ...
            'coordinate', dim(id).coord(:), ...
            'source', 'free axis');
    else
        info.indices.(field) = struct( ...
            'index', fixedIndex(id), ...
            'coordinate', fixedValue(id), ...
            'source', fixedSource{id});
    end
end
end

function dim = dimensionInfo(p)
dim = struct('id', {}, 'name', {}, 'coord', {}, 'aliases', {}, ...
    'defaultMode', {});

dim(1).id = 1;
dim(1).name = 'rho_x';
dim(1).coord = p.relative.rhoX.cells(:);
dim(1).aliases = {'rho_x', 'rhox', 'rhoX', 'xi_x', 'xrel'};
dim(1).defaultMode = 'nearest-zero';

dim(2).id = 2;
dim(2).name = 'rho_y';
dim(2).coord = p.relative.rhoY.cells(:);
dim(2).aliases = {'rho_y', 'rhoy', 'rhoY', 'xi_y', 'yrel'};
dim(2).defaultMode = 'nearest-zero';

dim(3).id = 3;
dim(3).name = 'X';
dim(3).coord = p.dg.X.nodes(:);
dim(3).aliases = {'X', 'x', 'chi'};
dim(3).defaultMode = 'middle-index';

dim(4).id = 4;
dim(4).name = 'Y';
dim(4).coord = p.dg.Y.nodes(:);
dim(4).aliases = {'Y', 'y'};
dim(4).defaultMode = 'middle-index';
end

function freeIds = parseFreeDimensions(freeDims, dim)
if ischar(freeDims) || isstring(freeDims)
    freeText = cellstr(freeDims);
    if isscalar(freeText) && contains(freeText{1}, ',')
        freeText = strtrim(strsplit(freeText{1}, ','));
    end
elseif iscell(freeDims)
    freeText = freeDims;
else
    freeText = [];
end

if isnumeric(freeDims)
    freeIds = freeDims(:).';
else
    freeIds = zeros(1, numel(freeText));
    for k = 1:numel(freeText)
        freeIds(k) = dimensionId(freeText{k}, dim);
    end
end

if numel(freeIds) ~= 2 || numel(unique(freeIds)) ~= 2 ...
        || any(freeIds < 1) || any(freeIds > 4)
    error('DG:Full2D:InvalidFreeDims', ...
        'freeDims must select exactly two different axes from X, Y, rho_x and rho_y.');
end
end

function id = dimensionId(name, dim)
name = char(name);
for k = 1:numel(dim)
    candidates = [{dim(k).name}, dim(k).aliases];
    if any(strcmpi(name, candidates))
        id = dim(k).id;
        return
    end
end

error('DG:Full2D:UnknownSliceDimension', ...
    'Unknown dimension "%s". Use X, Y, rho_x or rho_y.', name);
end

function [idx, value, source] = resolveFixedIndex(fixed, dim)
indexField = findMatchingField(fixed, indexCandidates(dim));
coordField = findMatchingField(fixed, coordinateCandidates(dim));

if ~isempty(indexField)
    idx = fixed.(indexField);
    validateattributes(idx, {'numeric'}, ...
        {'scalar', 'integer', 'positive', '<=', numel(dim.coord)}, ...
        mfilename, indexField);
    source = ['index field ', indexField];
elseif ~isempty(coordField)
    coordValue = fixed.(coordField);
    validateattributes(coordValue, {'numeric'}, {'scalar', 'finite'}, ...
        mfilename, coordField);
    [~, idx] = min(abs(dim.coord(:) - coordValue));
    source = ['nearest coordinate field ', coordField];
else
    switch dim.defaultMode
        case 'nearest-zero'
            [~, idx] = min(abs(dim.coord(:)));
            source = 'default nearest zero relative coordinate';
        otherwise
            idx = ceil(numel(dim.coord)/2);
            source = 'default middle center-coordinate index';
    end
end

value = dim.coord(idx);
end

function candidates = indexCandidates(dim)
aliases = [{dim.name}, dim.aliases];
candidates = cell(1, 2*numel(aliases));
for k = 1:numel(aliases)
    candidates{2*k-1} = [aliases{k}, 'Index'];
    candidates{2*k} = [aliases{k}, '_index'];
end
end

function candidates = coordinateCandidates(dim)
candidates = [{dim.name}, dim.aliases];
end

function field = findMatchingField(s, candidates)
field = '';
if ~isstruct(s)
    return
end

names = fieldnames(s);
for k = 1:numel(candidates)
    hit = find(strcmpi(names, candidates{k}), 1);
    if ~isempty(hit)
        field = names{hit};
        return
    end
end
end
