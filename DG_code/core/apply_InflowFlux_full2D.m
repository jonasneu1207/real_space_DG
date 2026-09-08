function flux = apply_InflowFlux_full2D(sideBoundary, rhoInside)
%APPLY_INFLOWFLUX_FULL2D Apply an already split Source/Drain boundary flux.
%
% With outward normal n, the characteristic upwind boundary flux is
%   F_n^* = A_n^+ rho_inside + A_n^- rho_boundary.
% Therefore outgoing components are never overwritten by reservoir data.

rhoInside = expandInsideState(rhoInside, size(sideBoundary.rhoBoundary));
rhoBoundary = effectiveBoundaryState(sideBoundary, rhoInside);
flux = sideBoundary.normalFlux.Aplus*rhoInside ...
     + sideBoundary.normalFlux.Aminus*rhoBoundary;
end

function rhoInside = expandInsideState(rhoInside, targetSize)
if isvector(rhoInside)
    rhoInside = rhoInside(:);
    if numel(rhoInside) ~= targetSize(1)
        error('DG:Full2D:InvalidInsideState', ...
            'Inside vector has %d entries, expected %d.', numel(rhoInside), targetSize(1));
    end
    rhoInside = repmat(rhoInside, 1, targetSize(2));
elseif ~isequal(size(rhoInside), targetSize)
    error('DG:Full2D:InvalidInsideState', ...
        'Inside state has size [%d,%d], expected [%d,%d].', ...
        size(rhoInside, 1), size(rhoInside, 2), targetSize(1), targetSize(2));
end
end

function rhoBoundary = effectiveBoundaryState(sideBoundary, rhoInside)
rhoBoundary = sideBoundary.rhoBoundary;
if ~isfield(sideBoundary, 'contactMask')
    return
end

contactMask = logical(sideBoundary.contactMask(:)).';
if numel(contactMask) ~= size(rhoBoundary, 2)
    error('DG:Full2D:InvalidContactMask', ...
        'Contact mask has %d entries, expected %d boundary columns.', ...
        numel(contactMask), size(rhoBoundary, 2));
end

closedFace = ~contactMask;
if ~any(closedFace)
    return
end

if isfield(sideBoundary, 'nonContactGhostOperator') ...
        && ~isempty(sideBoundary.nonContactGhostOperator)
    rhoBoundary(:, closedFace) = ...
        sideBoundary.nonContactGhostOperator*rhoInside(:, closedFace);
else
    rhoBoundary(:, closedFace) = rhoInside(:, closedFace);
end
end
