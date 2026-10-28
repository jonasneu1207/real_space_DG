function result = solve_self_consistent_qtbm2D(mat,poisson,options)
%SOLVE_SELF_CONSISTENT_QTBM2D Quantum transport with shared 2D Poisson.
% Convergence tests the actual input density/residual and UNMIXED update.
% Final phi and transport always refer to the same evaluated potential.
ramp = repmat(bias_ramp_qtbm2D(mat),1,mat.Ny);
phi = options.bias(1)*ramp; points = cell(1,numel(options.bias));
if ~isempty(options.initialPhi)
    validateattributes(options.initialPhi,{'double'},{'real','finite','size',[mat.Nx,mat.Ny]});
    phi = options.initialPhi;
end
status = 'self-consistent-complete';
% Restriction of the SHARED Poisson operator to the Y-homogeneous subspace.
A = poisson.matrix(1:mat.Nx,:)*repmat(speye(mat.Nx),mat.Ny,1);
free = 2:mat.Nx-1; q = 1.602176634e-19;
chargeScale = q*max(mat.Nd(:));
if isfield(mat,'Na'), acceptors = mat.Na(:,1); else, acceptors = mat.N_a(:,1)*1e6; end
net = mat.Nd(:,1)-acceptors;
for k = 1:numel(options.bias)
    bias = options.bias(k);
    if k>1, phi = phi+(bias-options.bias(k-1))*ramp; end
    phi(1,:) = 0; phi(end,:) = bias;
    previousDensity = []; history = nan(options.maxIterations,5);
    converged = false; xHistory = []; rHistory = []; lastResidual = Inf;
    baseX = []; trialStep = []; backtracks = 0;
    inverseJacobian = []; previousX = []; previousF = [];
    for iteration = 1:options.maxIterations
        potentialEV = mat.Vc-phi;
        contacts = get_contacts_qtbm2D(mat,potentialEV,options);
        transport = solve_transport_qtbm2D(mat,potentialEV,contacts,options);
        [nPoisson,projection] = project_density_qtbm2D(transport,mat);
        state = struct('phi',phi,'bias',bias,'nPoisson',nPoisson);
        [candidate,pInfo] = poisson_step_qtbm2D(state,poisson);
        if ~pInfo.converged, error('QTBM2D:PoissonLinearSolve','Poisson linear update failed.'); end
        densityChange = Inf;
        if ~isempty(previousDensity)
            densityChange = max(abs(nPoisson-previousDensity),[],'all')/ ...
                max([max(abs(nPoisson(:))),max(mat.Nd(:)),1]);
        end
        history(iteration,:) = [pInfo.unmixedStep,densityChange, ...
            pInfo.physicalRelresBefore,transport.currentDensity, ...
            transport.diagnostics.energy.evaluations];
        if options.verbose && (iteration==1 || mod(iteration,10)==0)
            fprintf('QTBM V=%.4g, it=%d: dphi=%.3g V, dn=%.3g, Poisson=%.3g\n', ...
                bias,iteration,pInfo.unmixedStep,densityChange,pInfo.physicalRelresBefore);
        end
        % Backtrack extrapolations that worsen the fixed-point residual.
        % Rejected trials do not enter the secant history or density test.
        if ~isempty(baseX) && pInfo.unmixedStep>1.1*lastResidual ...
                && backtracks<8 && iteration<options.maxIterations
            trialStep = trialStep/2;
            phi = repmat(baseX+trialStep,1,mat.Ny);
            phi(1,:) = 0; phi(end,:) = bias;
            backtracks = backtracks+1;
            continue
        end
        converged = pInfo.unmixedStep<=options.potentialTolerance ...
            && densityChange<=options.densityTolerance ...
            && pInfo.physicalRelresBefore<=options.poissonResidualTolerance;
        if converged || iteration==options.maxIterations, break; end
        previousDensity = nPoisson;
        % Couple transport to the shared approximate-Newton Poisson map.
        % Broyden uses secants of the physical residual; Anderson uses the map.
        x = phi(:,1); r = mean(candidate-phi,2);
        f = (A(free,:)*x-q*(nPoisson(free,1)-net(free)))/chargeScale;
        if pInfo.unmixedStep>3*lastResidual
            xHistory = []; rHistory = [];
        end
        if strcmp(options.couplingMethod,'broyden')
            if isempty(inverseJacobian)
                J = A(free,free);
                if strcmp(options.poissonMethod,'newton-raphson')
                    J = J-spdiags(q^2*nPoisson(free,1)/(1.380649e-23*mat.Temp),0,numel(free),numel(free));
                end
                inverseJacobian = full(J\(chargeScale*speye(numel(free))));
            elseif ~isempty(previousX)
                s = x(free)-previousX; y = f-previousF;
                By = inverseJacobian*y; denominator = s.'*By;
                if abs(denominator)>1e-8*norm(s)*norm(By)
                    inverseJacobian = inverseJacobian+(s-By)*(s.'*inverseJacobian)/denominator;
                end
            end
            previousX = x(free); previousF = f;
            step = zeros(mat.Nx,1);
            step(free) = -options.mixing*inverseJacobian*f;
        else
            xHistory = [xHistory,x]; rHistory = [rHistory,r]; %#ok<AGROW>
            keep = max(1,size(xHistory,2)-options.andersonDepth):size(xHistory,2);
            xHistory = xHistory(:,keep); rHistory = rHistory(:,keep);
            step = options.mixing*r;
            if size(xHistory,2)>1
                dx = diff(xHistory,1,2); dr = diff(rHistory,1,2);
                gamma = pinv(dr,1e-6*max(norm(dr,2),realmin))*r;
                step = step-(dx+options.mixing*dr)*gamma;
            end
        end
        step = step*min(1,options.maxPotentialStepV/max(norm(step,Inf),realmin));
        baseX = x; trialStep = step; backtracks = 0;
        phi = repmat(x+step,1,mat.Ny);
        lastResidual = pInfo.unmixedStep;
        phi(1,:) = 0; phi(end,:) = bias;
    end
    points{k} = struct('bias',bias,'phi',phi,'potentialEV',potentialEV, ...
        'transport',transport,'nPoisson',nPoisson,'projection',projection, ...
        'converged',converged,'iterations',iteration,'poissonInfo',pInfo, ...
        'history',array2table(history(1:iteration,:),'VariableNames', ...
        {'unmixedStepV','densityChange','poissonResidual','currentDensity','energyEvaluations'}));
    if options.verbose
        fprintf('QTBM V=%g V: converged=%d, iterations=%d, j=%g A/m^2\n', ...
            bias,converged,iteration,transport.currentDensity);
    end
    if ~converged
        status = 'self-consistency-not-converged'; points = points(1:k);
        warning('QTBM2D:SelfConsistency','Stopped at V=%g V; invalid IV entry is NaN.',bias);
        break
    end
end
result = collect_sweep_qtbm2D(points,mat,options,status);
end
