function normalFlux = get_NormalFlux_full2D(p, normal)
%GET_NORMALFLUX_FULL2D Characteristic normal flux for the full 2D scaffold.
%
% For an outward normal n=(nX,nY), A_n = nX*A_X + nY*A_Y. Inflow means
% lambda(A_n) < 0 because the characteristic enters the physical X-Y domain
% against the outward normal. The boundary flux convention is
%
%   F_n^* = A_n^+ f_inside + A_n^- f_boundary.
%
% The rectangular Full-2D path only needs axis-aligned normals. For these
% faces the relative operators are separable, so the characteristic split can
% be formed from the one-dimensional rho_x or rho_y operator and lifted by a
% Kronecker product.

normal = normal(:).';
if numel(normal) ~= 2
    error('DG:Full2D:InvalidNormal', 'normal must contain [nX, nY].');
end

tol = 1e-12;
if isfield(p, 'boundary') && isfield(p.boundary, 'normalFluxTol')
    tol = p.boundary.normalFluxTol;
end

nX = normal(1);
nY = normal(2);
normalFlux = struct;
normalFlux.normal = normal;
normalFlux.signConvention = ...
    'A_n^+ multiplies the interior state; A_n^- multiplies prescribed inflow data.';

if abs(nY) <= tol
    [Aplus1D, Aminus1D, lambda1D] = splitCharacteristicMatrix(nX*p.relative.Ax1D, tol);
    normalFlux.An = nX*p.relative.Ax;
    normalFlux.Aplus = kron(speye(p.relative.NrhoY), Aplus1D);
    normalFlux.Aminus = kron(speye(p.relative.NrhoY), Aminus1D);
    normalFlux.lambda = repmat(lambda1D(:), p.relative.NrhoY, 1);
    normalFlux.axis = 'rho_x';
elseif abs(nX) <= tol
    [Aplus1D, Aminus1D, lambda1D] = splitCharacteristicMatrix(nY*p.relative.Ay1D, tol);
    normalFlux.An = nY*p.relative.Ay;
    normalFlux.Aplus = kron(Aplus1D, speye(p.relative.NrhoX));
    normalFlux.Aminus = kron(Aminus1D, speye(p.relative.NrhoX));
    normalFlux.lambda = kron(lambda1D(:), ones(p.relative.NrhoX, 1));
    normalFlux.axis = 'rho_y';
else
    maxDenseDof = 400;
    if isfield(p, 'boundary') && isfield(p.boundary, 'maxDenseNormalFluxDof')
        maxDenseDof = p.boundary.maxDenseNormalFluxDof;
    end
    if p.relative.nDof > maxDenseDof
        error('DG:Full2D:GeneralNormalTooLarge', ...
            ['General non-axis-aligned normal flux would require a dense ', ...
             '%d-by-%d characteristic split. Rectangular DG faces should use ', ...
             'axis-aligned normals.'], p.relative.nDof, p.relative.nDof);
    end
    normalFlux.An = nX*p.relative.Ax + nY*p.relative.Ay;
    [normalFlux.Aplus, normalFlux.Aminus, normalFlux.lambda] = ...
        splitCharacteristicMatrix(normalFlux.An, tol);
    normalFlux.axis = 'general';
end

normalFlux.inflowMask = real(normalFlux.lambda) < -tol;
normalFlux.outflowMask = real(normalFlux.lambda) > tol;
normalFlux.Aplus = sparse(normalFlux.Aplus);
normalFlux.Aminus = sparse(normalFlux.Aminus);
normalFlux.An = sparse(normalFlux.An);
end

function [Aplus, Aminus, lambda] = splitCharacteristicMatrix(A, tol)
[Phi, Lambda] = eig(full(A));
lambda = diag(Lambda);
lambdaPlus = lambda;
lambdaMinus = lambda;
lambdaPlus(real(lambdaPlus) <= tol) = 0;
lambdaMinus(real(lambdaMinus) >= -tol) = 0;
PhiInv = Phi\eye(size(Phi, 1));
Aplus = Phi*diag(lambdaPlus)*PhiInv;
Aminus = Phi*diag(lambdaMinus)*PhiInv;
Aplus(abs(Aplus) < 1e-12) = 0;
Aminus(abs(Aminus) < 1e-12) = 0;
end
