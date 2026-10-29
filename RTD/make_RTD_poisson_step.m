function [step,setup] = make_RTD_poisson_step(mat,method)
%MAKE_RTD_POISSON_STEP Connect material-grid Poisson to the DG outer loop.
% Builds the matrix once. Damping belongs to the OUTER mixing option only.
arguments
    mat (1,1) struct
    method {mustBeTextScalar} = 'direct'
end
if ~isfield(mat,'poisson'), mat.poisson = struct; end
if ~isfield(mat.poisson,'opt'), mat.poisson.opt = struct; end
mat.poisson.opt.solve = char(method);
mat.poisson.opt.alpha = 1;
if ~any(strcmpi(method,{'direct','newton-raphson','newton-rhapson'}))
    error('RTD:UnknownPoissonMethod','Use direct or newton-raphson.');
end
xBoundary = get_RTD_poisson_x_boundary(mat);
neumannX = strcmp(xBoundary,'neumann');
if neumannX && strcmpi(method,'direct')
    error('RTD:NeumannNeedsNewton','Neumann X requires newton-raphson.');
end
if isfield(mat,'boundary') && ~isempty(mat.boundary)
    expected = false(mat.Nx,mat.Ny);
    if ~neumannX, expected([1,end],:) = true; end
    if ~isequal(logical(mat.boundary),expected)
        error('RTD:PoissonContactModel','The coupled RTD supports full source/drain contacts only.');
    end
end
if isfield(mat,'Vbc') && ~isempty(mat.Vbc)
    error('RTD:PoissonContactModel','The coupled RTD sets contact values from bias; omit mat.Vbc.');
end
[M,ids,types,setup] = eval_Poisson2D_rtd(mat);
get_RTD_doping_SI(mat); % fail on inconsistent units before an expensive solve
periodic = strcmp(setup.yBoundary,'periodic');
q = mat.dg.params;
transportPeriodic = strcmp(get_YBoundaryType_full2D(mat),'periodic');
if periodic ~= transportPeriodic
    error('RTD:PoissonBoundaryMismatch','Poisson and transport must use the same Y periodicity.');
end
if isfield(q,'full2D_contactModeIntegrateKz') && ~q.full2D_contactModeIntegrateKz
    error('RTD:PoissonDensityUnits','Poisson needs volume density: enable full2D_contactModeIntegrateKz.');
end
setup.method = char(method);
setup.approximateNewton = ~strcmpi(method,'direct');
step = @poissonStep;
    function [phi,info] = poissonStep(state,~)
        [n,transfer] = project_DG_density_RTD(state.transport.nDG,state.transport.p,mat,periodic);
        % Never silently remove negative charge-density artifacts.
        % if min(n(:)) < 0
        %     error('RTD:NegativeElectronDensity', ...
        %         'Poisson received negative control-volume density (%g m^-3). Check the converged DG solution/grid.',min(n(:)));
        % end
        [phi,~,info] = solve_Poisson2D_rtd(M,types,ids,n,mat,state.phi,state.bias);
        info.chargeTransfer = transfer;
        info.electronDensity = n;
    end
end
