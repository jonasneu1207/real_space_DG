function [n, jx, jy, info] = get_Observables_full2D(rho, p, mat)
%GET_OBSERVABLES_FULL2D Observable scaffold for full physical 2D Wigner-DG.
%
% The final implementation should map rho(X,Y,eta_x,eta_y) to carrier
% density n(X,Y) and current densities jx(X,Y), jy(X,Y).

n = [];
jx = [];
jy = [];
info = struct;
info.currentComponents = {'jx', 'jy'};

error('DG:Full2D:NotImplemented', ...
    ['get_Observables_full2D is a placeholder for density and current ', ...
     'extraction on rho size [%d, %d] and grid [%d, %d], with %d parameter ', ...
     'fields and %d prepared outputs.'], ...
     size(rho, 1), size(rho, 2), mat.Nx, mat.Ny, numel(fieldnames(p)), ...
     numel(n) + numel(jx) + numel(jy) + numel(info.currentComponents));
end
