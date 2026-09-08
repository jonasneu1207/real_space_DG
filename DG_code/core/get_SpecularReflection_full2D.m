function reflection = get_SpecularReflection_full2D(p)
%GET_SPECULARREFLECTION_FULL2D Specular wall reflection in rho basis.
%
% Physical meaning: at Y-bottom/Y-top there is no reservoir. The normal
% transverse momentum is reflected, (k_x,k_y)->(k_x,-k_y), while k_x is kept.
% In the rho representation this is the discrete flip rho_y -> -rho_y:
%
%   rho_reflected(rho_x,rho_y) = rho_inside(rho_x,-rho_y).
%
% Masked non-contact pieces of the X-left/X-right outer faces are closed in
% the same sense. Their normal momentum is k_x, so the ghost state uses
% rho_x -> -rho_x while rho_y is kept.
%
% The relative vector order is (rho_x, rho_y), with rho_x fastest, so the
% sparse reflection matrices are
%   R_x = I_rho_y kron Flip_rho_x,
%   R_y = Flip_rho_y kron I_rho_x.

Ix = speye(p.relative.NrhoX);
Iy = speye(p.relative.NrhoY);
flipX = sparse(1:p.relative.NrhoX, p.relative.NrhoX:-1:1, 1, ...
    p.relative.NrhoX, p.relative.NrhoX);
flipY = sparse(1:p.relative.NrhoY, p.relative.NrhoY:-1:1, 1, ...
    p.relative.NrhoY, p.relative.NrhoY);
Rx = kron(Iy, flipX);
Ry = kron(flipY, Ix);

reflection = struct;
reflection.type = 'specular-reflection';
reflection.axis = 'rho_y';
reflection.R_x = Rx;
reflection.R_y = Ry;
reflection.x.axis = 'rho_x';
reflection.x.ghostOperator = Rx;
reflection.x.note = 'Vertical closed walls reflect k_x, represented as rho_x -> -rho_x.';
reflection.y.axis = 'rho_y';
reflection.y.ghostOperator = Ry;
reflection.y.note = 'Horizontal closed walls reflect k_y, represented as rho_y -> -rho_y.';
reflection.bottom.face = 'Y-bottom';
reflection.bottom.normal = p.domain.normals.YBottom;
reflection.bottom.ghostOperator = Ry;
reflection.top.face = 'Y-top';
reflection.top.normal = p.domain.normals.YTop;
reflection.top.ghostOperator = Ry;
reflection.note = ...
    'Y boundaries use reflected ghost data, not Dirichlet-zero and not reservoir injection.';
end
