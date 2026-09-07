function flux = apply_InflowFlux_full2D(sideBoundary, rhoInside)
%APPLY_INFLOWFLUX_FULL2D Apply an already split Source/Drain boundary flux.
%
% With outward normal n, the characteristic upwind boundary flux is
%   F_n^* = A_n^+ rho_inside + A_n^- rho_boundary.
% Therefore outgoing components are never overwritten by reservoir data.

rhoInside = expandInsideState(rhoInside, size(sideBoundary.rhoBoundary));
flux = sideBoundary.normalFlux.Aplus*rhoInside ...
     + sideBoundary.normalFlux.Aminus*sideBoundary.rhoBoundary;
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
