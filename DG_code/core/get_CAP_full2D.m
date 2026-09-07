function cap = get_CAP_full2D(p)
%GET_CAP_FULL2D Separable CAP on the two relative-coordinate axes.
%
% The CAP is placed only at the outer rho_x/rho_y boundaries, not on the
% physical X/Y device boundaries. Its separable discrete form is
%
%   C_rho = C_x kron I_y + I_x kron C_y.
%
% With MATLAB vector order (rho_x fastest), this is implemented as
%   C_rho = kron(I_y, C_x) + kron(C_y, I_x).

Wx = get_CAP(p.relative.rhoX.cells, p.relative.rhoX.length, p.relative.rhoX.delta);
Wy = get_CAP(p.relative.rhoY.cells, p.relative.rhoY.length, p.relative.rhoY.delta);
Wx = Wx(:);
Wy = Wy(:);

Cx = spdiags(Wx, 0, p.relative.NrhoX, p.relative.NrhoX);
Cy = spdiags(Wy, 0, p.relative.NrhoY, p.relative.NrhoY);
Ix = speye(p.relative.NrhoX);
Iy = speye(p.relative.NrhoY);
Crho = kron(Iy, Cx) + kron(Cy, Ix);

cap = struct;
cap.Cx = Cx;
cap.Cy = Cy;
cap.Crho = sparse(Crho);
cap.potentialMatrix = -1i*cap.Crho;
cap.profileX = Wx;
cap.profileY = Wy;
cap.profile2D = full(reshape(diag(cap.Crho), p.relative.NrhoX, p.relative.NrhoY));
cap.maskX = Wx ~= 0;
cap.maskY = Wy ~= 0;
cap.mask2D = cap.profile2D ~= 0;
cap.note = ...
    'CAP is a relative-coordinate absorber and must be assembled separately from physical boundary fluxes.';
end
