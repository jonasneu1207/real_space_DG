function wall = get_HardWallGeometry_full2D(p)
%GET_HARDWALLGEOMETRY_FULL2D Nodal Dirichlet domain in (rho_y,Y).
%
% A quantum wall constrains BOTH particle coordinates y1=Y+rho_y/2 and
% y2=Y-rho_y/2. Keep only samples strictly inside the resulting diamond.
% This is a nodal embedded-domain approximation on the existing DG/FV grid,
% NOT a boundary-fitted cut-cell discretization: between-node wall geometry
% requires a joint Y/rho_y refinement study. No extra absorbing penalty is
% used to implement the wall.
%
% testProjection = E*(E'*M*E)^(-1)*E'*M, with E injecting active unknowns.
% The exact (non-lumped) DG mass is used. Simply masking the strong operator
% on both sides would not preserve its weak/Galerkin structure at cut DOFs.

rhoY = p.relative.rhoY.cells(:);
Y = p.dg.Y.nodes(:).';
scale = max(abs([p.domain.Y(:);rhoY]));
tolerance = 64*eps(scale);
distance = min(Y-p.domain.Y(1),p.domain.Y(2)-Y);
active = abs(rhoY)/2 < distance-tolerance;
if ~any(active,'all')
    error('DG:Full2D:HardWallEmptyDomain', ...
        'No Y/rho_y samples lie inside the hard walls. Refine these grids.');
end

nR = numel(rhoY);
nY = numel(Y);
nLocal = p.dg.Y.nLocal;
M = full(p.dg.Y.M);
T = spalloc(nR*nY,nR*nY,nR*nY*nLocal);
for element = 1:p.dg.Y.nElements
    yIds = (element-1)*nLocal+(1:nLocal);
    for ir = 1:nR
        keep = active(ir,yIds);
        if any(keep)
            ids = ir+nR*(yIds-1);
            T(ids(keep),ids) = M(keep,keep)\M(keep,:);
        end
    end
end

wall = struct;
wall.type = 'two-point-dirichlet-hard-wall';
wall.boundsY = p.domain.Y;
wall.active = active;
wall.testProjection = T;
wall.arraySize = p.index.arraySize;
wall.activeDof = nnz(active)*p.relative.NrhoX*p.dg.X.nDof;
wall.constrainedDof = p.index.nTotal-wall.activeDof;
wall.tolerance = tolerance;
wall.discretization = 'nodal embedded diamond with exact-DG-mass restriction';
wall.note = ['Zero on/outside either particle-coordinate wall at all ', ...
    'stored samples. The oblique wall between samples is not fitted; ', ...
    'this is not a claim of exact continuum hard-wall dynamics.'];
end
