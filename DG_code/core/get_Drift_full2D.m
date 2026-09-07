function [G, info] = get_Drift_full2D(mat, p, Vxy)
%GET_DRIFT_FULL2D Potential-term scaffold for full physical 2D Wigner-DG.
%
% The nonlocal Wigner potential will later use
%   V(X+rho_x/2,Y+rho_y/2) - V(X-rho_x/2,Y-rho_y/2).
% For this first boundary-focused step we only provide the relative CAP
% contribution. The CAP is not mixed with Source/Drain or specular-reflection
% boundary fluxes.

G = [];
cap = get_CAP_full2D(p);

info = struct;
info.fullPotentialAssembled = false;
info.deviceSize = size(Vxy);
info.materialType = mat.type;
info.relativeCoordinates = p.relativeCoordinates;
info.CAP = cap;
info.CAPMatrixSize = size(cap.Crho);
info.CAPNnz = nnz(cap.Crho);
end
