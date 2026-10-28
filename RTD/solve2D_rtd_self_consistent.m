function results = solve2D_rtd_self_consistent(solver,mat,options)
%SOLVE2D_RTD_SELF_CONSISTENT Warm-started DG transport / electrostatic Poisson.
% Built-in Poisson: options.poissonMethod='direct' or 'newton-raphson'.
% Newton uses a LOCAL BOLTZMANN APPROXIMATION, not an exact quantum Jacobian.
% Optionally override the built-in step with options.poissonStep:
%   [phiCandidate,info] = poissonStep(state,mat)
% phiCandidate: real finite [mat.Nx,mat.Ny] electrostatic potential in VOLTS.
% info.converged: scalar logical, convergence of that Poisson step.
% state.phi, state.bias, state.iteration, state.biasIndex describe the input.
% state.transport is the full DG result. state.n is sampled on
% state.densityX [m] x state.densityY [m], NOT necessarily on mat.x x mat.y.
% Built-in Poisson conservatively transfers raw transport.nDG in both axes.
%
% Options: maxIterations=30, mixing=0.1, potentialTolerance=1e-5 V,
% densityTolerance=1e-4, poissonResidualTolerance=1e-5 (both relative),
% initialPhi=[], initialRho=[]. Built-in convergence also requires the actual
% Poisson residual at the CURRENT potential and quantum-transport density.
% Convergence uses the UNMIXED potential defect and density change. Returned
% rho, phi and Vc always belong to the same last transport evaluation.
% Only converged transport states are handed to the callback or reused.
% A converged previous bias also seeds the next bias on the same fixed grid.
if nargin < 3, options = struct; end
options = checkedOptions(options);
if ~strcmpi(solver,'DG')
    error('RTD:UnsupportedSolver','The self-consistency scaffold uses the full2D DG solver.');
end
rtdDir = fileparts(mfilename('fullpath'));
originalPath = path;
restorePath = onCleanup(@() path(originalPath));
addpath(fullfile(fileparts(rtdDir),'DG_code'),fullfile(fileparts(rtdDir),'DG_code','core'));
if ~isfield(mat,'dg') || ~isfield(mat.dg,'params')
    error('RTD:MissingDGParameters','Configure mat.dg.params before starting the loop.');
end
if isfield(mat.dg.params,'full2D_solver') && ~strcmpi(mat.dg.params.full2D_solver,'gmres')
    error('RTD:SelfConsistencyNeedsGmres','Set full2D_solver=''gmres'' to use iteration warm starts.');
end
mat.dg.params.full2D_solver = 'gmres';
mat.dg.params.full2D_solve = true;
validateattributes(mat.Vc,{'double','single'},{'real','finite','size',[mat.Nx,mat.Ny]});
validateattributes(mat.Vdd,{'double','single'},{'real','finite','nonempty','vector'});
bias = mat.Vdd(:).';
x = mat.x(:);
rampLength = x(end)-mat.Lc1-mat.Lc2;
if rampLength <= 0
    error('RTD:ContactGeometry','Contact regions must leave a positive central length.');
end
ramp = repmat(min(1,max(0,(x-mat.Lc1)/rampLength)),1,mat.Ny);
periodicY = strcmp(get_YBoundaryType_full2D(mat),'periodic');
useBuiltInPoisson = isempty(options.poissonStep);
poissonSetup = struct('method','external-callback');
if useBuiltInPoisson
    method = options.poissonMethod;
    if isempty(method)
        method = 'direct';
        if isfield(mat,'poisson') && isfield(mat.poisson,'opt') && isfield(mat.poisson.opt,'solve')
            method = mat.poisson.opt.solve;
        end
    end
    %% make Poisson step
    [options.poissonStep,poissonSetup] = make_RTD_poisson_step(mat,method);
    %%
end
rhoGuess = options.initialRho;
previousPhi = [];
poissonInfo = struct;
histories = cell(1,numel(bias));
sweep = struct('requestedBias',bias);
for ib = 1:numel(bias)
    if ib == 1
        phi = bias(ib)*ramp;
        if ~isempty(options.initialPhi), phi = options.initialPhi; end
    else
        phi = previousPhi+(bias(ib)-bias(ib-1))*ramp;
    end
    checkPotential(phi,mat,bias(ib),periodicY);
    previousDensity = [];
    history = repmat(struct('iteration',0,'transportFlag',NaN, ...
        'transportRelres',NaN,'warmStartUsed',false,'initialTrueRelres',NaN, ...
        'gmresIterations',0,'densityChange',Inf,'potentialDefect',Inf, ...
        'poissonConverged',false,'poissonPhysicalRelres',NaN),options.maxIterations,1);
    status = 'iteration-limit';
    converged = false;
    for iteration = 1:options.maxIterations
        poissonInfo = struct; % no stale Poisson diagnostics after transport failure
        V = mat.Vc-phi; % electron energy in eV, electrostatic phi in volts
        efL = mat.Vc(1,1)+0.05;
        efR = mat.Vc(end,1)+0.05-bias(ib);
        if isfield(mat,'EfL'), efL = mat.EfL; end
        if isfield(mat,'EfR'), efR = mat.EfR-bias(ib); end
        %% solve transport
        transport = solve_transport_DG_full2D(mat,V,efL,efR,rhoGuess);
        %%
        if isempty(transport.rho)
            error('RTD:SolveSkipped','Transport returned no state; self-consistency cannot continue.');
        end
        h = history(iteration);
        h.iteration = iteration;
        h.transportFlag = transport.info.solve.flag;
        h.transportRelres = transport.info.solve.relres;
        h.warmStartUsed = transport.info.solve.initialGuess.used;
        h.initialTrueRelres = transport.info.solve.initialGuess.initialTrueRelres;
        h.gmresIterations = transport.info.solve.progress.completedIterations;
        history(iteration) = h;
        if h.transportFlag ~= 0 || ~isfinite(h.transportRelres) ...
                || h.transportRelres > transport.info.solve.tolerance
            status = 'transport-not-converged';
            break
        end
        if ~isempty(previousDensity)
            h.densityChange = norm(transport.n(:)-previousDensity(:)) ...
                /max([norm(transport.n(:)),norm(previousDensity(:)),realmin]);
        end
        state = struct('phi',phi,'bias',bias(ib),'biasIndex',ib, ...
            'iteration',iteration,'transport',transport,'n',transport.n, ...
            'densityX',transport.chi,'densityY',transport.y, ...
            'materialX',mat.x(:)*transport.p.coordinateScale, ...
            'materialY',mat.y(:)*transport.p.coordinateScale);
        %% Poisson lösen 
        % Built-in step or explicitly supplied external callback.
        [candidate,poissonInfo] = options.poissonStep(state,mat);
        %%
        checkPotential(candidate,mat,bias(ib),periodicY);
        if ~isstruct(poissonInfo) || ~isscalar(poissonInfo) ...
                || ~isfield(poissonInfo,'converged') ...
                || ~islogical(poissonInfo.converged) || ~isscalar(poissonInfo.converged)
            error('RTD:InvalidPoissonInfo','poissonStep must return info.converged as a scalar logical.');
        end
        h.poissonConverged = poissonInfo.converged;
        h.potentialDefect = max(abs(candidate(:)-phi(:)));
        poissonResidualOK = true;
        if useBuiltInPoisson
            h.poissonPhysicalRelres = poissonInfo.physicalRelresBefore;
            poissonResidualOK = isfinite(h.poissonPhysicalRelres) ...
                && h.poissonPhysicalRelres <= options.poissonResidualTolerance;
        end
        history(iteration) = h;
        if ~poissonInfo.converged
            status = 'poisson-not-converged';
            break
        end
        if h.potentialDefect <= options.potentialTolerance ...
                && h.densityChange <= options.densityTolerance && poissonResidualOK
            status = 'fixed-point-converged';
            converged = true;
            break
        end
        if iteration == options.maxIterations
            break % do NOT return an updated potential with a stale transport state
        end
        previousDensity = transport.n;
        rhoGuess = transport.rho; % full PHYSICAL solution, no preconditioned coordinates
        phi = phi+options.mixing*(candidate-phi);
    end
    histories{ib} = history(1:iteration);
    sweep.bias(ib) = bias(ib);
    sweep.current(:,ib) = transport.j;
    sweep.density(:,:,ib) = transport.n;
    sweep.potential(:,:,ib) = V;
    sweep.phi(:,:,ib) = phi;
    sweep.flag(ib) = transport.info.solve.flag;
    sweep.relres(ib) = transport.info.solve.relres;
    sweep.converged(ib) = converged;
    sweep.iterations(ib) = iteration;
    sweep.status{ib} = status;
    if ~converged
        break % no bias continuation with an unconverged state
    end
    previousPhi = phi;
    rhoGuess = transport.rho;
end
results = transport;
results.transportStatus = transport.status;
results.status = ['self-consistency-',status];
results.phi = phi;
results.sweep = sweep;
results.Id = mean(sweep.current,1);
results.n_electron = sweep.density;
results.Vc = sweep.potential;
results.mat = mat;
results.poisson = poissonInfo;
results.selfConsistency = struct('converged',all(sweep.converged) && ib==numel(bias), ...
    'status',status,'history',{histories(1:ib)},'iterations',sweep.iterations, ...
    'poissonImplementedHere',useBuiltInPoisson,'poissonSetup',poissonSetup, ...
    'poissonCallback',func2str(options.poissonStep), ...
    'note','Built-in: fully ionized electron-only Poisson; approximate Newton uses local Boltzmann screening. External callbacks remain caller responsibility.');
end

function options = checkedOptions(input)
if ~isstruct(input) || ~isscalar(input)
    error('RTD:InvalidSelfConsistencyOptions','options must be a scalar struct.');
end
options = struct('poissonStep',[],'poissonMethod','','maxIterations',100,'mixing',0.1, ...
    'potentialTolerance',5e-4,'densityTolerance',5e-4,'poissonResidualTolerance',8e-2, ...
    'initialPhi',[],'initialRho',[]);                                                            % Default Werte für die Poisson-Loop
names = fieldnames(input);                                                                       % Namen alle Felder die durch Anwender definiert wurden in main
for i = 1:numel(names)                                                                           % Überschreibung der Default-Wete durch Angaben vom Anwender
    if ~isfield(options,names{i})
        error('RTD:UnknownSelfConsistencyOption','Unknown option: %s.',names{i});
    end
    options.(names{i}) = input.(names{i});                                                       %
end
if ~isempty(options.poissonStep) && ~isa(options.poissonStep,'function_handle')
    error('RTD:InvalidPoissonCallback','poissonStep must be empty (built-in) or a function handle.');
end
validateattributes(options.poissonResidualTolerance,{'numeric'},{'real','finite','scalar','positive'});
validateattributes(options.maxIterations,{'numeric'},{'real','finite','scalar','integer','positive'});
validateattributes(options.mixing,{'numeric'},{'real','finite','scalar','>',0,'<=',1});
validateattributes(options.potentialTolerance,{'numeric'},{'real','finite','scalar','positive'});
validateattributes(options.densityTolerance,{'numeric'},{'real','finite','scalar','positive'});
end

function checkPotential(phi,mat,bias,periodicY)
if ~isfloat(phi) || ~isreal(phi) || ~isequal(size(phi),[mat.Nx,mat.Ny]) ...
        || any(~isfinite(phi(:))) || isa(phi,'gpuArray')
    error('RTD:InvalidPoissonPotential','phi must be a finite real CPU array [mat.Nx,mat.Ny] in volts.');
end
tolerance = 1e-10*max(1,abs(bias));
if max(abs(phi(1,:))) > tolerance || max(abs(phi(end,:)-bias)) > tolerance
    error('RTD:PoissonContactPotential','This RTD scaffold requires phi(source)=0 and phi(drain)=bias.');
end
if periodicY && max(abs(phi(:,1)-phi(:,end))) > tolerance
    error('RTD:NonperiodicPoissonPotential','The first and last material Y columns of phi must agree.');
end
end

