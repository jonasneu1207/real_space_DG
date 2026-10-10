function pre = get_LineTransversePreconditioner_full2D(mat, p, Vxy, sys, useGPU)
%GET_LINETRANSVERSEPRECONDITIONER Fixed approximate-factorization inverse.
%   PRE.apply(r) approximates A\r without assembling A. Only this auxiliary
%   operator is transformed to eigenmodes of -i*D_rhoX. For each mode it
%   solves a shifted upwind X line and a coupled (rhoY,Y) problem:
%       M = (sigma*I + Xref)*(sigma*I + Tref)/sigma.
%   Tref retains the frozen-X transverse BDD operator, FV potential and CAP.
%   Xref uses the scalar reference mass. X-dependent coefficients, X BDD
%   corrections, closed contact-face reflection and off-mode couplings are
%   left to the OUTER Krylov iteration. No physical coefficient is modified.
%
%   Setup is CPU double. CPU applications use sparse LU factors; GPU
%   applications use dense LU pages and pagemldivide, with no host transfers.
%   The factor-memory guard is conservative for CPU and checked before setup.
%   See full2D_line_transverse.md for controls, costs and limitations.

params = mat.dg.params;
storageDefault = 'sparse';
if useGPU
    storageDefault = 'dense-pages';
end
storage = char(option(params,'full2D_lineFactorStorage',storageDefault));
if ~any(strcmp(storage,{'sparse','dense-pages'})) || (useGPU && strcmp(storage,'sparse'))
    error('DG:Full2D:LineStorage', ...
        'Use sparse or dense-pages on CPU; GPU requires dense-pages.');
end
usePages = strcmp(storage,'dense-pages');
sz = p.index.arraySize;
nModes = sz(1);
nT = sz(2)*sz(4);
nLocal = p.dg.X.nLocal;
nElem = p.dg.X.nElements;
factorBytes = 32*nT^2*nModes + 32*nLocal^2*nModes*nElem ...
    + 16*nModes^2 + 8*nT*nModes;
maxBytes = option(params, 'full2D_lineMaxFactorBytes', 4*1024^3);
if ~isscalar(maxBytes) || ~isfinite(maxBytes) || maxBytes <= 0
    error('DG:Full2D:LineMemoryLimit', 'lineMaxFactorBytes must be positive and finite.');
end
if factorBytes > maxBytes
    error('DG:Full2D:LineMemoryLimit', ...
        'Estimated factors %.3g GiB exceed full2D_lineMaxFactorBytes %.3g GiB.', ...
        factorBytes/1024^3, maxBytes/1024^3);
end
referenceIndex = option(params, 'full2D_lineReferenceXIndex', 1);
if ~isscalar(referenceIndex) || referenceIndex < 1 ...
        || referenceIndex > mat.Nx || referenceIndex ~= fix(referenceIndex)
    error('DG:Full2D:LineReferenceIndex', 'Reference X index must be an integer in 1:mat.Nx.');
end
setupTimer = tic;
[U, lambda] = eig(full(p.relative.Ax1D), 'vector');
lambda = real(lambda);
lambda(abs(lambda) < 100*eps(max(abs(lambda)))) = 0;
[K, V] = transverseReference(mat, p, Vxy, sys, referenceIndex);
capY = repmat(sys.drift.CAP.profileY(:), sz(4), 1) ...
    * sys.drift.scale.capCoefficient;
capX = real(sum(conj(U).*(sys.drift.CAP.profileX(:).*U), 1)) ...
    * sys.drift.scale.capCoefficient;
if strcmp(sys.drift.potentialDiscretization, 'fv-consistent')
    P = spdiags(ones(nModes+1,1)*[1,1], [-1,0], nModes+1, nModes);
    potentialWeight = real(sum(conj(U).*((P'*P/4)*U),1));
else
    potentialWeight = ones(1,nModes);
end
Qx = sys.diff.transportScales.Qx;
xRate = abs(Qx)*max(abs(lambda))/min(p.dg.X.jacobian);
tRate = norm(K,inf)+norm(V,inf)+max(abs(capY))+max(abs(capX));
relativeShift = option(params, 'full2D_lineRelativeShift', 0.003);
sigma = option(params, 'full2D_lineShift', relativeShift*max([xRate,tRate,realmin]));
if ~isscalar(sigma) || ~isreal(sigma) || ~isfinite(sigma) || sigma <= 0
    error('DG:Full2D:LineShift', 'The preconditioner shift must be positive, real and finite.');
end
data = struct('sz',sz, 'U',U, 'nLocal',nLocal, 'nElem',nElem, ...
    'sigma',sigma, 'usePages',usePages);
[data.xInverse, data.xCoupling, data.forward] = xFactors(p.dg.X, Qx*lambda, sigma);
wall = [];
if isfield(sys,'hardWall')
    wall = sys.hardWall;
end
data.L = cell(nModes,1);
data.R = cell(nModes,1);
data.permutation = zeros(nT,nModes);
for k = 1:nModes
    B = K + potentialWeight(k)*V + spdiags(capY+capX(k),0,nT,nT);
    if ~isempty(wall)
        active = spdiags(double(wall.active(:)),0,nT,nT);
        B = wall.testProjection*B*active ...
            + wall.constraintScale*(speye(nT)-active);
    end
    B = B/sigma + speye(nT);
    if usePages
        % Dense page solves need triangular factors and the row permutation.
        [Ld,Rd,permutation] = lu(full(B),'vector');
        data.L{k} = Ld;
        data.R{k} = Rd;
        data.permutation(:,k) = permutation(:);
    else
        data.L{k} = decomposition(B,'lu');
    end
end
if usePages
    data.L = cat(3,data.L{:});
    data.R = cat(3,data.R{:});
end
if useGPU
    data.L = gpuArray(data.L);
    data.R = gpuArray(data.R);
    data.permutation = gpuArray(data.permutation);
    data.U = gpuArray(data.U);
    data.xInverse = gpuArray(data.xInverse);
    data.xCoupling = gpuArray(data.xCoupling);
    wait(gpuDevice);
end
pre = struct('enabled',true, 'apply',@(r) applyInverse(data,r));
pre.info = struct('requested','line-transverse', 'method','line-transverse', ...
    'enabled',true, 'nTotal',p.index.nTotal, 'shift',sigma, ...
    'referenceXIndex',referenceIndex, 'transverseDof',nT, 'modeCount',nModes, ...
    'factorBytesUpperBound',factorBytes, 'factorStorage',storage, ...
    'scratchBytesEstimate',8*16*p.index.nTotal, ...
    'setupSeconds',toc(setupTimer), 'executionDevice','cpu', ...
    'fixedLinear',true, 'precision','double', ...
    'note',['Frozen-X transverse BDD/potential/CAP and reference-mass upwind X lines; ', ...
    'X BDD corrections and off-characteristic couplings remain in the outer operator.']);
if useGPU
    pre.info.executionDevice = 'gpu';
end
end

function [K,V] = transverseReference(mat,p,Vxy,sys,ix)
% A two-rhoX-cell, two-X-element strip is sufficient: Qx=0 and all
% coefficients are X invariant. Reuse the actual quadrature/BDD builders.
ref = mat;
ref.V = repmat(reshape(Vxy(ix,:),1,[]),mat.Nx,1);
massNames = {'me_x','me_y'};
for j = 1:numel(massNames)
    name = massNames{j};
    if isfield(ref,name) && ~isscalar(ref.(name))
        field = reshape(ref.(name),mat.Nx,mat.Ny);
        ref.(name) = reshape(repmat(field(ix,:),mat.Nx,1),size(ref.(name)));
    end
end
q = ref.dg.params;
q.N_rho_x = 2;
q.N_K_X = 2;
% Existing DG builders require at least two elements for their face stencil.
q.N_el_X = 2;
q.full2D_gpu = false;
q.full2D_Q_diff_x = 0;
q.full2D_cap_scale = 0;
q.full2D_assembleDiffMatrix = true;
q.full2D_assembleDriftMatrix = true;
q.full2D_reservoirModel = 'zero';
% No reservoir construction is necessary for a homogeneous operator block.
q.full2D_sourceRho = zeros(2*p.relative.NrhoY,p.dg.Y.nDof);
q.full2D_drainRho = q.full2D_sourceRho;
ref.dg.params = q;
pr = initParams_full2D(ref,ref.V,[],[]);
boundary = get_Boundary_full2D(ref,pr,[],[],ref.V);
if strcmp(sys.massModel,'position-dependent-bdd')
    Kfull = get_Diff_variableMass_full2D(ref,pr,boundary,false);
else
    Kfull = get_Diff_full2D(ref,pr,boundary,false);
end
Vfull = get_Drift_full2D(ref,pr,ref.V);
ids = reshape(1:pr.index.nTotal,pr.index.arraySize);
ids = reshape(ids(1,:,1,:),[],1);
K = Kfull(ids,ids);
V = Vfull(ids,ids);
if strcmp(sys.drift.potentialDiscretization,'fv-consistent')
    % The extracted rhoX diagonal of P_x'*P_x/4 equals 1/2.
    V = 2*V;
end
end

function [inverse,coupling,forward] = xFactors(axis,speeds,sigma)
n = axis.nLocal;
Pleft = zeros(n); Pleft(1,1) = 1;
Pright = zeros(n); Pright(end,end) = 1;
PfromLeft = zeros(n); PfromLeft(1,end) = 1;
PfromRight = zeros(n); PfromRight(end,1) = 1;
K1 = full(axis.M\(axis.S+0.5*Pleft-0.5*Pright));
K2 = full(axis.M\(0.5*Pleft+0.5*Pright));
inverse = complex(zeros(n,n,numel(speeds),axis.nElements));
coupling = inverse;
forward = speeds >= 0;
for e = 1:axis.nElements
    for k = 1:numel(speeds)
        a = speeds(k)/(sigma*axis.jacobian(e));
        diagonal = eye(n) + a*K1+abs(a)*K2;
        inverse(:,:,k,e) = diagonal\eye(n);
        if forward(k)
            coupling(:,:,k,e) = -a*full(axis.M\PfromLeft);
        else
            coupling(:,:,k,e) = a*full(axis.M\PfromRight);
        end
    end
end
end

function z = applyInverse(d,r)
% Transform only rhoX; arrange [X-local, transverse, mode, X-element].
sz = d.sz;
w = d.U'*reshape(r,sz(1),[]);
w = permute(reshape(w,sz),[3,2,4,1]);
w = reshape(w,d.nLocal,d.nElem,sz(2)*sz(4),sz(1));
w = permute(w,[1,3,4,2]);
v = zeros(size(w),'like',w);
groups = {find(d.forward),find(~d.forward)};
orders = {1:d.nElem,d.nElem:-1:1};
for g = 1:2
    modes = groups{g};
    if isempty(modes)
        continue
    end
    previous = zeros(d.nLocal,sz(2)*sz(4),numel(modes),'like',w);
    for e = orders{g}
        block = w(:,:,modes,e) ...
            - pagemtimes(d.xCoupling(:,:,modes,e),previous);
        previous = pagemtimes(d.xInverse(:,:,modes,e),block);
        v(:,:,modes,e) = previous;
    end
end
% [transverse, X-DG, mode], with factors reused for every X right-hand side.
v = permute(v,[2,1,4,3]);
v = reshape(v,sz(2)*sz(4),sz(3),sz(1));
if d.usePages
    permuted = zeros(size(v),'like',v);
    for k = 1:sz(1)
        permuted(:,:,k) = v(d.permutation(:,k),:,k);
    end
    v = pagemldivide(d.R,pagemldivide(d.L,permuted));
else
    for k = 1:sz(1)
        v(:,:,k) = d.L{k}\v(:,:,k);
    end
end
v = reshape(v,sz(2),sz(4),sz(3),sz(1));
v = permute(v,[4,1,3,2]);
z = reshape(d.U*reshape(v,sz(1),[]),[],1)/d.sigma;
end

function value = option(params,name,default)
value = default;
if isfield(params,name) && ~isempty(params.(name))
    value = params.(name);
end
end
