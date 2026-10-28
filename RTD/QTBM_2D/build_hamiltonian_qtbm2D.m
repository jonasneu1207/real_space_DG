function hamiltonian = build_hamiltonian_qtbm2D(mat,potentialEV)
%BUILD_HAMILTONIAN_QTBM2D Constant-mass longitudinal FD Hamiltonian [eV].
% H = Hx + Ey + Ez; do not include transverse energies in Hx again.
x = double(mat.x(:))*1e-9; dx = diff(x);
if numel(x)<3 || any(dx<=0) || max(abs(dx-dx(1)))>1e-10*dx(1)
    error('QTBM2D:Grid','At least three uniformly spaced X nodes are required.');
end
if ~isequal(size(potentialEV),[mat.Nx,mat.Ny]) || any(~isfinite(potentialEV(:))) ...
        || ~isreal(potentialEV) || max(abs(potentialEV-potentialEV(:,1)),[],'all')>1e-10
    error('QTBM2D:NonseparablePotential','Potential must be finite, real and Y homogeneous.');
end
m = mat.me_x;
if any(m(:)<=0) || any(~isfinite(m(:))) || max(abs(m(:)-m(1)))>1e-12*m(1)
    error('QTBM2D:Mass','Only globally constant positive longitudinal mass is supported.');
end
hbar = 1.054571817e-34; q = 1.602176634e-19; m0 = 9.1093837015e-31;
t = hbar^2/(2*m(1)*m0*dx(1)^2*q);
v = double(potentialEV(:,1)); N = numel(x);
H = spdiags([-t*ones(N,1),v+2*t,-t*ones(N,1)],-1:1,N,N);
hamiltonian = struct('matrixEV',H,'potentialEV',v,'hoppingEV',t, ...
    'diagonalEV',v+2*t,'xM',x,'dxM',dx(1),'size',N, ...
    'sourceEdgeEV',v(1),'drainEdgeEV',v(end));
end
