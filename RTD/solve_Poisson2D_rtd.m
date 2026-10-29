function [phi,iter_err,info] = solve_Poisson2D_rtd(M,node_type,node_num,n,mat,phi0,bias)
%SOLVE_POISSON2D_RTD Electrostatic direct or approximate Newton update [V].
% M approximates +div(epsilon*grad(phi)); hence M*phi=q*(n-Nd+Na).
% n [m^-3] is a nodal-control-volume average on the material grid.
% Nd/Na are SI; N_d/N_a are legacy cm^-3 (checked if both are present).
% mat.poisson.opt.solve: 'direct' (default), 'newton-raphson', or legacy
% 'newton-rhapson'. Newton is ONE approximate outer Newton step with
% dn/dphi=q*n/(kB*T), NOT the full quantum-transport Jacobian.
% opt.alpha damps the update (default 1); the coupled runner mixes externally.
% info.converged means the LINEAR UPDATE was solved, not self-consistency.
% physicalRelresBefore tests Poisson with the actual input transport density.
% Pure Neumann X uses Newton screening with fixed reservoir chemical
% potentials. No gauge pin and no subtraction of net charge are applied.
% Direct fixed-charge inversion is deliberately rejected for this mode.
arguments
    M {mustBeFloat,mustBeReal,mustBeFinite}
    node_type {mustBeNumeric,mustBeReal,mustBeFinite}
    node_num {mustBeNumeric,mustBeReal,mustBeFinite}
    n {mustBeFloat,mustBeReal,mustBeFinite}                             %{mustBeFloat,mustBeReal,mustBeFinite,mustBeNonnegative}
    mat (1,1) struct
    phi0 {mustBeFloat,mustBeReal,mustBeFinite}
    bias (1,1) {mustBeFloat,mustBeReal,mustBeFinite}
end
q = 1.602176634e-19;
kB = 1.380649e-23;
Nx = mat.Nx; Ny = mat.Ny; N = Nx*Ny;
validateattributes(M,{'double','single'},{'size',[N,N]});
validateattributes(n,{'double','single'},{'size',[Nx,Ny]});
validateattributes(phi0,{'double','single'},{'size',[Nx,Ny]});
[Nd,Na] = get_RTD_doping_SI(mat);
xBoundary = get_RTD_poisson_x_boundary(mat);
neumannX = strcmp(xBoundary,'neumann');
if neumannX && isfield(mat,'Vbc') && ~isempty(mat.Vbc)
    error('RTD:PoissonBoundaryConflict','Do not supply fixed mat.Vbc with Neumann X.');
end
if isequal(size(node_num),[Nx+2,Ny+2]) && isequal(size(node_type),[Nx+2,Ny+2])
    ids = node_num(2:Nx+1,2:Ny+1);
    types = node_type(2:Nx+1,2:Ny+1);
else
    ids = node_num; types = node_type;
end
if ~isequal(ids,reshape(1:N,Nx,Ny)) || ~isequal(size(types),[Nx,Ny])
    error('RTD:PoissonNodeMap','Rebuild the node maps with eval_Poisson2D_rtd.');
end
periodic = any(types(:)==3);
Nu = Ny-double(periodic);
if periodic && (~all(types(:,end)==3) || any(any(types(:,1:Nu)==3)))
    error('RTD:PoissonNodeMap','Only the final Y column may contain periodic duplicates.');
end
bc = zeros(Nx,Ny); bc(end,:) = bias;
if isfield(mat,'Vbc') && ~isempty(mat.Vbc), bc = mat.Vbc; end
validateattributes(bc,{'double','single'},{'real','finite','size',[Nx,Ny]});
if periodic
    requirePeriodic(n,'electron density');
    requirePeriodic(Nd,'donor density'); requirePeriodic(Na,'acceptor density');
    requirePeriodic(bc,'contact potential'); requirePeriodic(phi0,'initial potential');
end
uniqueTypes = types(:,1:Nu);
free = find(uniqueTypes(:)==1); fixed = find(uniqueTypes(:)==2);
if numel(free)+numel(fixed) ~= Nx*Nu || (~neumannX && isempty(fixed)) ...
        || (neumannX && ~isempty(fixed))
    error('RTD:PoissonNodeMap','Invalid or ungrounded Poisson node types.');
end
% Reduce the duplicate Y endpoint, avoiding an extra physical cell/charge.
map = reshape(1:Nx*Nu,Nx,Nu);
if periodic, map = [map,map(:,1)]; end
expand = sparse(1:N,map(:),1,N,Nx*Nu);
A = double(M(1:Nx*Nu,:))*expand;
old = double(phi0(:,1:Nu)); old = old(:);
bc = double(bc(:,1:Nu)); bc = bc(:);
nf = double(n(:,1:Nu)); nf = nf(:);
net = Nd(:,1:Nu)-Na(:,1:Nu); net = net(:);
charge = q*(nf-net);
v = old; v(fixed) = bc(fixed);
Av = A*v;
F = Av-charge;
scale = max([norm(q*nf(free),Inf),norm(q*net(free),Inf),norm(Av(free),Inf),q]);
mode = 'direct'; alpha = 1;
if isfield(mat,'poisson') && isfield(mat.poisson,'opt')
    opt = mat.poisson.opt;
    if isfield(opt,'solve'), mode = lower(char(opt.solve)); end
    if isfield(opt,'alpha'), alpha = opt.alpha; end
end
validateattributes(alpha,{'numeric'},{'real','finite','scalar','>',0,'<=',1});
Aff = A(free,free);
switch mode
    case 'direct'
        if neumannX
            error('RTD:NeumannNeedsNewton', ...
                'Neumann X requires newton-raphson and fixed reservoir Fermi levels, not direct fixed-charge inversion.');
        end
        delta = -(Aff\F(free));
        J = Aff;
        derivativeModel = 'fixed transport density';
    case {'newton-raphson','newton-rhapson'}
        validateattributes(mat.Temp,{'numeric'},{'real','finite','scalar','positive'});
        if neumannX && (any(nf<0) || ~any(nf>0))
            error('RTD:NeumannDensityResponse', ...
                'Neumann Newton needs nonnegative electron density with a positive screening response. Check the transport discretization.');
        end
        J = Aff-spdiags(q^2*nf(free)/(kB*mat.Temp),0,numel(free),numel(free));
        delta = -(J\F(free));
        derivativeModel = 'local Boltzmann approximate quantum-density Jacobian';
    otherwise
        error('RTD:UnknownPoissonMethod','Use direct or newton-raphson, not %s.',mode);
end
linearRelres = norm(J*delta+F(free),Inf)/max(norm(F(free),Inf),q);
target = v; target(free) = v(free)+delta;
updated = v; updated(free) = v(free)+alpha*delta;
phi = reshape(expand*updated,Nx,Ny);
iter_err = norm(updated(free)-old(free));
after = A*updated-charge;
% Physical charge per unit unresolved z length [C/m]. Count the periodic
% endpoint once and use HALF X control volumes, including for diagnostics.
wx = ones(Nx,1)*mat.dx*1e-9; wx([1,end]) = wx([1,end])/2;
wy = ones(1,Nu)*mat.dy*1e-9;
if ~periodic && Nu>1, wy([1,end]) = wy([1,end])/2; end
volumes = wx*wy; volumes = volumes(:);
netCharge = -sum(volumes.*charge);
chargeMagnitude = max(q*max(sum(volumes.*abs(net)),sum(volumes.*abs(nf))),realmin);
info = struct('converged',all(isfinite(updated)) && isfinite(linearRelres) && linearRelres<1e-8, ...
    'xBoundary',xBoundary,'neumannCompatibilityRequired',neumannX, ...
    'netChargePerLength',netCharge,'relativeChargeImbalance',abs(netCharge)/chargeMagnitude, ...
    'method',mode,'derivativeModel',derivativeModel,'linearRelres',linearRelres, ...
    'physicalRelresBefore',norm(F(free),Inf)/scale, ...
    'physicalRelresAfter',norm(after(free),Inf)/scale, ...
    'physicalResidualBefore',norm(F(free),Inf),'residualScale',scale, ...
    'unmixedStep',norm(delta,Inf),'alpha',alpha, ...
    'candidateUnmixed',reshape(expand*target,Nx,Ny), ...
    'note','converged refers to the linear update; recompute quantum density for outer convergence.');
end

function requirePeriodic(value,label)
if max(abs(value(:,1)-value(:,end))) > 1e-10*max([1;abs(value(:))])
    error('RTD:NonperiodicPoissonData','The %s must agree at the included Y endpoints.',label);
end
end

