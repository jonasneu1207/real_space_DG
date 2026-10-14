function defect = get_PeriodicYStateError_full2D(u, p)
%GET_PERIODICYSTATEERROR_FULL2D Two-particle periodicity on the double cover.
% Accept a full state or a contact face [rhoX*rhoY,nY]. No projection is
% performed: this diagnostic exposes nonphysical half-period components.
if numel(u) == p.relative.nDof*p.dg.Y.nDof
    values = reshape(u,p.relative.NrhoX,p.relative.NrhoY,1,p.dg.Y.nDof);
else
    values = reshape(u,p.index.arraySize);
end
shifted = values(:,p.periodicY.relativeShift,:,p.periodicY.centerShift);
defect = norm(values(:)-shifted(:))/max(norm(values(:)),realmin);
end
