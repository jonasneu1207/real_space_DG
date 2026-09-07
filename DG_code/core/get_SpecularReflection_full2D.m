function reflection = get_SpecularReflection_full2D(p)
%GET_SPECULARREFLECTION_FULL2D Specular Y-boundary reflection in rho basis.
%
% Physical meaning: at Y-bottom/Y-top there is no reservoir. The normal
% transverse momentum is reflected, (k_x,k_y)->(k_x,-k_y), while k_x is kept.
% In the rho representation this is the discrete flip rho_y -> -rho_y:
%
%   rho_reflected(rho_x,rho_y) = rho_inside(rho_x,-rho_y).
%
% The relative vector order is (rho_x, rho_y), with rho_x fastest, so the
% sparse reflection matrix is R_y = Flip_rho_y kron I_rho_x.

Ix = speye(p.relative.NrhoX);
flipY = sparse(1:p.relative.NrhoY, p.relative.NrhoY:-1:1, 1, ...
    p.relative.NrhoY, p.relative.NrhoY);
Ry = kron(flipY, Ix);

reflection = struct;
reflection.type = 'specular-reflection';
reflection.axis = 'rho_y';
reflection.R_y = Ry;
reflection.bottom.face = 'Y-bottom';
reflection.bottom.normal = p.domain.normals.YBottom;
reflection.bottom.ghostOperator = Ry;
reflection.top.face = 'Y-top';
reflection.top.normal = p.domain.normals.YTop;
reflection.top.ghostOperator = Ry;
reflection.note = ...
    'Y boundaries use reflected ghost data, not Dirichlet-zero and not reservoir injection.';
end
