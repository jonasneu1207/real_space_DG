function [inflow, info] = get_InflowBoundary_full2D(mat, p, EfL, EfR, Vxy)
%GET_INFLOWBOUNDARY_FULL2D Source/drain inflow data for full 2D Wigner-DG.
%
% X-left is the Source face with outward normal (-1,0). X-right is the Drain
% face with outward normal (1,0). For both faces only A_n^- is prescribed by
% reservoir data; A_n^+ remains the interior outflow contribution.

if nargin < 5
    Vxy = [];
end

sourceFlux = get_NormalFlux_full2D(p, p.domain.normals.XLeft);
drainFlux = get_NormalFlux_full2D(p, p.domain.normals.XRight);

[sourceRho, sourceInfo] = getBoundaryRhoData(mat, p, 'source', p.dg.Y.nDof, EfL, Vxy);
[drainRho, drainInfo] = getBoundaryRhoData(mat, p, 'drain', p.dg.Y.nDof, EfR, Vxy);

inflow = struct;
inflow.source = makeInflowSide('X-left', 'Source', p.domain.normals.XLeft, ...
    sourceFlux, sourceRho, p.dg.faces.global.left.centerDofs, EfL, sourceInfo);
inflow.drain = makeInflowSide('X-right', 'Drain', p.domain.normals.XRight, ...
    drainFlux, drainRho, p.dg.faces.global.right.centerDofs, EfR, drainInfo);

info = struct;
info.source = sourceInfo;
info.drain = drainInfo;
info.signConvention = sourceFlux.signConvention;
end

function side = makeInflowSide(face, reservoir, normal, normalFlux, rhoBoundary, centerDofs, Ef, dataInfo)
side = struct;
side.type = 'characteristic-inflow';
side.face = face;
side.reservoir = reservoir;
side.normal = normal;
side.fermiLevel = Ef;
side.normalFlux = normalFlux;
side.rhoBoundary = rhoBoundary;
side.centerDofs = centerDofs(:);
side.dataInfo = dataInfo;
side.inflowComponents = nnz(normalFlux.inflowMask);
side.outflowComponents = nnz(normalFlux.outflowMask);
end

function [rhoBoundary, dataInfo] = getBoundaryRhoData(mat, p, sideName, nFaceDof, Ef, Vxy)
params = struct;
if isfield(mat, 'dg') && isfield(mat.dg, 'params')
    params = mat.dg.params;
end

directName = ['full2D_', sideName, 'Rho'];
longitudinalName = ['full2D_', sideName, 'RhoX'];
phaseName = ['full2D_', sideName, 'PhaseX'];

dataInfo = struct;
dataInfo.side = sideName;
dataInfo.faceDof = nFaceDof;

if isfield(params, directName) && ~isempty(params.(directName))
    rhoRel = params.(directName);
    dataInfo.origin = directName;
elseif isfield(params, longitudinalName) && ~isempty(params.(longitudinalName))
    opts = struct;
    opts.inputBasis = 'rho_x';
    opts.rhoYProfile = readOptionalVector(params, 'full2D_rhoYProfile');
    [rhoRel, transformInfo] = transform_InflowToRho2D_full2D(params.(longitudinalName), p, opts);
    dataInfo.origin = longitudinalName;
    dataInfo.transform = transformInfo;
elseif isfield(params, phaseName) && ~isempty(params.(phaseName))
    opts = struct;
    opts.inputBasis = 'phase_x';
    opts.rhoYProfile = readOptionalVector(params, 'full2D_rhoYProfile');
    [rhoRel, transformInfo] = transform_InflowToRho2D_full2D(params.(phaseName), p, opts);
    dataInfo.origin = phaseName;
    dataInfo.transform = transformInfo;
else
    [rhoRel, reservoirInfo] = get_DefaultReservoirRho_full2D(mat, p, sideName, Ef, Vxy);
    dataInfo.origin = reservoirInfo.origin;
    dataInfo.defaultReservoir = reservoirInfo;
end

rhoBoundary = expandFaceData(rhoRel, p.relative.nDof, nFaceDof);
dataInfo.size = size(rhoBoundary);
end

function data = readOptionalVector(params, name)
data = [];
if isfield(params, name)
    data = params.(name);
end
end

function rhoBoundary = expandFaceData(rhoData, nRelative, nFaceDof)
if isvector(rhoData)
    rhoData = rhoData(:);
    if numel(rhoData) ~= nRelative
        error('DG:Full2D:InvalidBoundaryData', ...
            'Boundary vector has %d entries, expected %d relative DOFs.', ...
            numel(rhoData), nRelative);
    end
    rhoBoundary = repmat(rhoData, 1, nFaceDof);
else
    if ~isequal(size(rhoData), [nRelative, nFaceDof])
        error('DG:Full2D:InvalidBoundaryData', ...
            'Boundary matrix has size [%d,%d], expected [%d,%d].', ...
            size(rhoData, 1), size(rhoData, 2), nRelative, nFaceDof);
    end
    rhoBoundary = rhoData;
end
end
