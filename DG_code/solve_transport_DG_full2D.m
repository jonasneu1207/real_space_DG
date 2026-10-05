function DG = solve_transport_DG_full2D(mat, Vxy, EfL, EfR)
%SOLVE_TRANSPORT_DG_FULL2D Stationary Full-2D Wigner-DG/FV rho solver.
%
% Coordinates:
%   X, Y          center-of-mass coordinates, rectangular DG
%   rho_x, rho_y relative coordinates, finite volumes
%
% Global unknown:
%   F(iRhoX, iRhoY, iX, iY), vectorized as F(:). rho_x is fastest, then
%   rho_y, X and Y.
%
% This routine now builds the Full-2D stationary operator in rho basis:
%   - sparse/Kronecker DG transport and boundary fluxes,
%   - FV-consistent potential blocks in rho_x/rho_y and separable CAP on
%     the outer relative-coordinate boundaries.
%
% Matrix policy:
%   Small systems are assembled as sparse matrices and solved directly by
%   default. Large systems keep A empty and expose sysInfo.apply(u). A
%   matrix-free Krylov solve is available via params.full2D_solve = true
%   and params.full2D_solver = 'bicgstab' or 'gmres'.
%
% Krylov progress:
%   params.full2D_showSolverProgress = true enables a progress display for
%   BICGSTAB and GMRES. params.full2D_solverProgressMode may be 'auto'
%   (GUI waitbar in the MATLAB desktop, text otherwise), 'waitbar', 'text'
%   or 'off'. The display is throttled by
%   params.full2D_solverProgressUpdateSeconds (default 0.25 s).
%
% GPU execution:
%   params.full2D_gpu = true runs BICGSTAB or GMRES with gpuArray data and
%   a fully device-resident matrix-free operator. No automatic GPU check or
%   device selection is performed. The selected/default MATLAB GPU is used,
%   and the converged rho vector is gathered before observables are built.

solverDir = fileparts(mfilename('fullpath'));
addpath(fullfile(solverDir, 'core'));

if nargin < 4
    EfR = [];
end
if nargin < 3
    EfL = [];
end
if nargin < 2 || isempty(Vxy)
    Vxy = mat.V;
end

p = initParams_full2D(mat, Vxy, EfL, EfR);
[A, rhs, sysInfo] = get_SysM_full2D(mat, p, Vxy, EfL, EfR);
[rho, solveInfo, A, rhs] = solveStationarySystem(A, rhs, sysInfo, mat, p);

DG = struct;
DG.status = solveInfo.status;
DG.basis = 'rho-full2D';
DG.flux = sysInfo.diff.fluxType;
DG.fluxX = sysInfo.diff.fluxTypeX;
DG.fluxY = sysInfo.diff.fluxTypeY;
DG.yBoundaryType = sysInfo.diff.yBoundaryType;
DG.coordinates = p.index.order;
DG.centerDiscretization = p.centerDiscretization;
DG.relativeDiscretization = p.relativeDiscretization;
DG.p = p;
DG.A = A;
DG.rhs = rhs;
DG.rho = rho;
DG.boundary = sysInfo.boundaryData;
DG.info = sysInfo;
DG.info.solve = solveInfo;
DG.chi = p.dg.X.nodes(:);
DG.yDG = p.dg.Y.nodes(:);
DG.rho_x = p.relative.rhoX.cells(:);
DG.rho_y = p.relative.rhoY.cells(:);
DG.xi = p.relative.rhoX.cells(:);
DG.getRhoSlice2D = @(freeDims, varargin) ...
    get_RhoSlice2D_full2D(rho, p, freeDims, varargin{:});
DG.slice2D = DG.getRhoSlice2D;
DG.info.slice2D = struct( ...
    'function', 'get_RhoSlice2D_full2D', ...
    'usage', ['[S,ax,si] = DG.getRhoSlice2D({''X'',''rho_x''}, ', ...
        'struct(''Y'',0,''rho_y'',0));'], ...
    'note', ['Slices are extracted from DG.rho using the global order ', ...
        'rho(iRhoX,iRhoY,iX,iY), so DG.RHO does not have to be stored.']);

if isempty(rho)
    DG.n = [];
    DG.j = [];
    DG.jx = [];
    DG.jy = [];
    DG.RHO = [];
    return
end

[n, jx, jy, observableInfo] = get_Observables_full2D(rho, p, mat);
DG.n = n;
DG.jx = jx;
DG.jy = jy;
DG.j = observableInfo.jxIntegratedOverY;
DG.nDG = observableInfo.nDG;
DG.jxDG = observableInfo.jxDG;
DG.jyDG = observableInfo.jyDG;
DG.y = observableInfo.y;
DG.info.observables = observableInfo;

if shouldStoreRHOArray(mat, p.index.nTotal)
    DG.RHO = reshape(rho, p.index.arraySize);
else
    DG.RHO = [];
    DG.info.solve.RHONote = ...
        'RHO array was not duplicated; use DG.rho and reshape with DG.p.index.arraySize, ordered as rho_x,rho_y,X,Y, if needed.';
end
end

function [rho, solveInfo, A, rhs] = solveStationarySystem(A, rhs, sysInfo, mat, p)
params = getDGParams(mat);
nTotal = p.index.nTotal;
mode = readParam(params, 'full2D_solve', 'auto');
maxAutoSolveDof = readParam(params, 'full2D_maxAutoSolveDof', 1e7);
maxMatrixFreeSolveDof = readParam(params, 'full2D_maxMatrixFreeSolveDof', 2e7);
useGPU = readLogicalParam(params, 'full2D_gpu', false);

rho = [];
solveInfo = struct;
solveInfo.requestedMode = mode;
solveInfo.nTotal = nTotal;
solveInfo.matrixAssembled = sysInfo.fullMatrixAssembled;
solveInfo.matrixFreeAvailable = sysInfo.matrixFreeAvailable;
solveInfo.tolerance = readParam(params, 'full2D_solverTol', 1e-4);
solveInfo.maxIterations = readParam(params, 'full2D_solverMaxIt', 1400);

[doSolve, skipStatus] = shouldSolve(mode, nTotal, maxAutoSolveDof);
if ~doSolve
    solveInfo.status = skipStatus;
    solveInfo.reason = sprintf(['Operator is ready but no solve was run. ', ...
        'nTotal=%d, full2D_solve=%s, full2D_maxAutoSolveDof=%d.'], ...
        nTotal, modeToString(mode), maxAutoSolveDof);
    return
end

if isempty(A) && nTotal > maxMatrixFreeSolveDof
    solveInfo.status = 'operator-ready-solve-skipped-matrix-free-limit';
    solveInfo.reason = sprintf(['Matrix-free solve skipped for %d DOFs. ', ...
        'Increase full2D_maxMatrixFreeSolveDof only when memory/runtime are acceptable.'], ...
        nTotal);
    return
end

solverName = chooseSolver(params, isempty(A), nTotal, useGPU);
solveInfo.solver = solverName;
solveInfo.gpu = struct( ...
    'requested', useGPU, ...
    'enabled', false, ...
    'resultGathered', false, ...
    'note', 'GPU setup has not been required for this solver path.');
if useGPU && strcmp(solverName, 'direct')
    error('DG:Full2D:GPUDirectSolverUnsupported', ...
        ['full2D_gpu=true supports full2D_solver=''bicgstab'' or ', ...
         '''gmres''. Select one of these matrix-free Krylov solvers.']);
end
preconditioner = buildPreconditioner(sysInfo, params, solverName, ...
    nTotal, useGPU);
solveInfo.preconditioner = preconditioner.info;

if strcmp(solverName, 'direct') && isempty(A)
    maxDirectSolveDof = readParam(params, 'full2D_maxDirectSolveDof', 150000);
    if nTotal > maxDirectSolveDof
        solveInfo.status = 'operator-ready-direct-solve-skipped-size';
        solveInfo.reason = sprintf(['Direct solve would require assembling ', ...
            '%d DOFs. Increase full2D_maxDirectSolveDof only for small ', ...
            'diagnostic systems, or use bicgstab/gmres.'], nTotal);
        return
    end
end

if isempty(rhs)
    rhs = sysInfo.assembleRhs();
end

applyA = [];
rhsSolve = rhs;
if any(strcmp(solverName, {'bicgstab', 'gmres'}))
    [applyA, rhsSolve, solveInfo.gpu] = prepareKrylovExecution( ...
        sysInfo, rhs, useGPU);
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch solverName
    case 'direct'
        if isempty(A)
            A = sysInfo.assemble();
        end
        rho = A\rhs;
        solveInfo.status = 'solved-direct';
        solveInfo.flag = 0;
        solveInfo.relres = norm(A*rho - rhs)/max(norm(rhs), eps);
        solveInfo.iterations = 1;

    case 'bicgstab'
        [rho, flag, relres, iter, resvec, progressInfo] = ...
            runBicgstabWithProgress(applyA, rhsSolve, preconditioner, ...
            solveInfo.tolerance, solveInfo.maxIterations, params, '');
        if shouldRetryBlockJacobiWithRowAbs(flag, preconditioner, params)
            firstAttempt = struct( ...
                'method', preconditioner.info.method, ...
                'flag', flag, ...
                'relres', relres, ...
                'iterations', iter, ...
                'progress', progressInfo);
            preconditioner = enableRowAbsPreconditioner(preconditioner, ...
                sysInfo, params, nTotal, ...
                sprintf(['BICGSTAB with Block-Jacobi returned flag %d; ', ...
                'retrying with row-absolute-sum scaling.'], flag));
            if preconditioner.enabled
                [rho, flag, relres, iter, resvec, progressInfo] = ...
                    runBicgstabWithProgress(applyA, rhsSolve, preconditioner, ...
                    solveInfo.tolerance, solveInfo.maxIterations, params, ...
                    'rowabs retry');
                preconditioner.info.retryFrom = firstAttempt;
                progressInfo.previousAttempt = firstAttempt.progress;
            end
        end
        solveInfo.status = iterativeStatus('bicgstab', flag);
        solveInfo.flag = flag;
        solveInfo.relres = relres;
        solveInfo.iterations = iter;
        solveInfo.resvec = resvec;
        solveInfo.progress = progressInfo;

    case 'gmres'
        if useGPU
            restartDefault = readParam(params, ...
                'full2D_gpuGmresRestart', 20);
        else
            restartDefault = [];
        end
        restart = readParam(params, 'full2D_gmresRestart', restartDefault);
        solveInfo.gpu.gmresRestart = restart;
        [rho, flag, relres, iter, resvec, progressInfo] = ...
            runGmresWithProgress(applyA, rhsSolve, preconditioner, restart, ...
            solveInfo.tolerance, solveInfo.maxIterations, nTotal, params);
        solveInfo.status = iterativeStatus('gmres', flag);
        solveInfo.flag = flag;
        solveInfo.relres = relres;
        solveInfo.iterations = iter;
        solveInfo.resvec = resvec;
        solveInfo.progress = progressInfo;

    otherwise
        error('DG:Full2D:UnknownSolver', ...
            'Unknown full2D_solver "%s". Use auto, direct, bicgstab or gmres.', ...
            solverName);
end
if useGPU && ~isempty(rho)
    rho = gather(rho);
    solveInfo.gpu.resultGathered = true;
end
solveInfo.preconditioner = preconditioner.info;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
function [applyA, rhsSolve, info] = prepareKrylovExecution(sysInfo, rhs, useGPU)
%PREPAREKRYLOVEXECUTION Place all iteration-resident data on CPU or GPU.
if ~useGPU
    applyA = @(u) sysInfo.apply(u);
    rhsSolve = rhs;
    info = struct( ...
        'requested', false, ...
        'enabled', false, ...
        'resultGathered', false, ...
        'matrixFreeOperator', true, ...
        'note', 'Krylov operator and vectors remain on the CPU.');
    return
end

if ~isfield(sysInfo, 'makeGPUApply') || isempty(sysInfo.makeGPUApply)
    error('DG:Full2D:MissingGPUOperator', ...
        'The Full-2D system does not expose makeGPUApply().');
end

% This intentionally performs no canUseGPU/gpuDevice query. gpuArray uses
% the device selected by the user (or MATLAB's current default) and reports
% the native MATLAB error if GPU support is unavailable.
applyA = sysInfo.makeGPUApply();
rhsSolve = gpuArray(rhs);
info = struct( ...
    'requested', true, ...
    'enabled', true, ...
    'resultGathered', false, ...
    'matrixFreeOperator', true, ...
    'underlyingClass', classUnderlying(rhsSolve), ...
    'note', ['RHS, Krylov vectors, transport factors, drift diagonal and ', ...
        'supported preconditioner data reside on the selected GPU.']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
function [x, flag, relres, iter, resvec, progressInfo] = ...
        runBicgstabWithProgress(applyA, rhs, preconditioner, tolerance, ...
        maxIterations, params, attemptLabel)
%RUNBICGSTABWITHPROGRESS Run BICGSTAB with a throttled progress callback.
%
% MATLAB's BICGSTAB has no public iteration callback. Its preconditioner is
% applied once per half iteration, so the monitored preconditioner provides
% the live iteration count without changing the Krylov recurrence. If no
% numerical preconditioner is requested, the monitored action is identity.

progress = createKrylovProgress(params, 'bicgstab', maxIterations, [], ...
    numel(rhs), attemptLabel);
cleanupProgress = onCleanup(@() progress.close());

if progress.enabled
    monitoredA = @(u) progress.applyOperator(applyA, u);
    monitoredM = @(r) progress.applyPreconditioner( ...
        preconditioner.apply, r);
    [x, flag, relres, iter, resvec] = bicgstab(monitoredA, rhs, ...
        tolerance, maxIterations, monitoredM);
elseif preconditioner.enabled
    [x, flag, relres, iter, resvec] = bicgstab(applyA, rhs, ...
        tolerance, maxIterations, preconditioner.apply);
else
    [x, flag, relres, iter, resvec] = bicgstab(applyA, rhs, ...
        tolerance, maxIterations);
end

flag = gatherIfGPU(flag);
relres = gatherIfGPU(relres);
iter = gatherIfGPU(iter);
resvec = gatherIfGPU(resvec);
progress.finish(iter, resvec, flag, relres);
progressInfo = progress.getInfo();
clear cleanupProgress
end

function [x, flag, relres, iter, resvec, progressInfo] = ...
        runGmresWithProgress(applyA, rhs, preconditioner, restart, ...
        tolerance, maxIterations, nTotal, params)
%RUNGMRESWITHPROGRESS Run GMRES and monitor inner Krylov iterations.
%
% With restart=[], MATLAB interprets maxIterations as the maximum number of
% inner iterations. With a finite restart, maxIterations counts outer
% restart cycles; the progress label therefore reports both the accumulated
% inner iterations and the current outer cycle.

progress = createKrylovProgress(params, 'gmres', maxIterations, restart, ...
    nTotal, '');
cleanupProgress = onCleanup(@() progress.close());

if progress.enabled
    monitoredA = @(u) progress.applyOperator(applyA, u);
    monitoredM = @(r) progress.applyPreconditioner( ...
        preconditioner.apply, r);
    [x, flag, relres, iter, resvec] = gmres(monitoredA, rhs, restart, ...
        tolerance, maxIterations, monitoredM);
elseif preconditioner.enabled
    [x, flag, relres, iter, resvec] = gmres(applyA, rhs, restart, ...
        tolerance, maxIterations, preconditioner.apply);
else
    [x, flag, relres, iter, resvec] = gmres(applyA, rhs, restart, ...
        tolerance, maxIterations);
end

flag = gatherIfGPU(flag);
relres = gatherIfGPU(relres);
iter = gatherIfGPU(iter);
resvec = gatherIfGPU(resvec);
progress.finish(iter, resvec, flag, relres);
progressInfo = progress.getInfo();
clear cleanupProgress
end

function progress = createKrylovProgress(params, solverName, ...
        maxIterations, restart, nTotal, attemptLabel)
%CREATEKRYLOVPROGRESS Create GUI/text progress callbacks for Krylov solvers.
%
% The live count is inferred from operator/preconditioner applications
% because MATLAB's public bicgstab/gmres interfaces expose no callback. The
% final display uses RESVEC to report all work actually performed. ITER is
% retained separately because MATLAB can return the index of an earlier,
% better approximation when a solver does not converge.

showProgress = readLogicalParam(params, ...
    'full2D_showSolverProgress', true);
requestedMode = lower(char(readParam(params, ...
    'full2D_solverProgressMode', 'auto')));
if ~showProgress || any(strcmp(requestedMode, ...
        {'off', 'none', 'false', 'no'}))
    mode = 'off';
elseif strcmp(requestedMode, 'auto')
    if usejava('desktop')
        mode = 'waitbar';
    else
        mode = 'text';
    end
elseif any(strcmp(requestedMode, {'waitbar', 'text'}))
    mode = requestedMode;
else
    error('DG:Full2D:UnknownSolverProgressMode', ...
        ['Unknown full2D_solverProgressMode "%s". Use auto, waitbar, ', ...
         'text or off.'], requestedMode);
end

updateSeconds = readParam(params, ...
    'full2D_solverProgressUpdateSeconds', 0.25);
updateSeconds = max(0, double(updateSeconds));
[maxWorkUnits, gmresInner, gmresOuter, gmresRestarted] = ...
    krylovWorkLimits(solverName, maxIterations, restart, nTotal);

preconditionerCalls = 0;
operatorCalls = 0;
updateCount = 0;
lastCompleted = 0;
lastUpdateTimer = tic;
waitbarHandle = [];
textLineOpen = false;
finalIterations = [];
finalFlag = [];
finalRelres = [];

displayName = upper(solverName);
if ~isempty(attemptLabel)
    displayName = sprintf('%s (%s)', displayName, attemptLabel);
end

if strcmp(mode, 'waitbar')
    try
        waitbarHandle = waitbar(0, progressMessage(0), ...
            'Name', 'Full-2D Wigner solver', 'NumberTitle', 'off');
    catch
        mode = 'text';
    end
end

progress = struct;
progress.enabled = ~strcmp(mode, 'off');
progress.applyOperator = @applyOperator;
progress.applyPreconditioner = @applyPreconditioner;
progress.finish = @finish;
progress.close = @closeProgress;
progress.getInfo = @getInfo;

    function y = applyOperator(baseOperator, u)
        y = baseOperator(u);
        operatorCalls = operatorCalls + 1;
        if strcmp(solverName, 'bicgstab') && preconditionerCalls > 0
            updateProgress(0.5*preconditionerCalls, false);
        end
    end

    function y = applyPreconditioner(basePreconditioner, r)
        if isempty(basePreconditioner)
            y = r;
        else
            y = basePreconditioner(r);
        end
        preconditionerCalls = preconditionerCalls + 1;
        if strcmp(solverName, 'gmres')
            completed = gmresLiveIterations(preconditionerCalls, ...
                gmresInner, gmresRestarted);
            updateProgress(completed, false);
        end
    end

    function updateProgress(completed, force)
        completed = min(max(double(completed), 0), maxWorkUnits);
        if ~force && updateCount > 0 && toc(lastUpdateTimer) < updateSeconds
            return
        end
        if ~force && completed <= lastCompleted
            return
        end

        lastCompleted = completed;
        lastUpdateTimer = tic;
        updateCount = updateCount + 1;
        fraction = min(completed/max(maxWorkUnits, 1), 1);
        message = progressMessage(completed);
        if strcmp(mode, 'waitbar')
            if isgraphics(waitbarHandle)
                waitbar(fraction, waitbarHandle, message);
                drawnow limitrate nocallbacks
            end
        elseif strcmp(mode, 'text')
            renderTextProgress(fraction, message, false);
        end
    end

    function finish(iter, resvec, flag, relres)
        finalIterations = iter;
        finalFlag = flag;
        finalRelres = relres;
        completed = exactCompletedIterations(solverName, resvec);
        lastCompleted = completed;
        updateCount = updateCount + 1;
        message = finalMessage(completed, flag, relres);
        if strcmp(mode, 'waitbar')
            if isgraphics(waitbarHandle)
                waitbar(1, waitbarHandle, message);
                drawnow limitrate nocallbacks
            end
        elseif strcmp(mode, 'text')
            renderTextProgress(1, message, true);
        end
    end

    function message = progressMessage(completed)
        if strcmp(solverName, 'gmres') && gmresRestarted
            cycle = min(max(ceil(max(completed, 1)/gmresInner), 1), ...
                gmresOuter);
            message = sprintf(['%s: %.0f / %.0f innere Iterationen ', ...
                '(Zyklus %d / %d)'], displayName, completed, ...
                maxWorkUnits, cycle, gmresOuter);
        else
            message = sprintf('%s: %g / %g Iterationen', displayName, ...
                completed, maxWorkUnits);
        end
    end

    function message = finalMessage(completed, flag, relres)
        if flag == 0
            state = 'konvergiert';
        else
            state = sprintf('beendet, flag=%d', flag);
        end
        message = sprintf('%s: %s nach %g Iterationen, relres=%.3e', ...
            displayName, state, completed, relres);
    end

    function renderTextProgress(fraction, message, finishLine)
        barWidth = 30;
        filled = min(floor(fraction*barWidth), barWidth);
        bar = [repmat('#', 1, filled), ...
            repmat('-', 1, barWidth-filled)];
        fprintf('\r[%s] %s                              ', bar, message);
        textLineOpen = true;
        if finishLine
            fprintf('\n');
            textLineOpen = false;
        end
    end

    function closeProgress()
        if isgraphics(waitbarHandle)
            delete(waitbarHandle);
        end
        if strcmp(mode, 'text') && textLineOpen
            fprintf('\n');
            textLineOpen = false;
        end
    end

    function info = getInfo()
        info = struct;
        info.enabled = ~strcmp(mode, 'off');
        info.mode = mode;
        info.solver = solverName;
        info.maxIterations = maxIterations;
        info.maxWorkUnits = maxWorkUnits;
        info.gmresRestart = restart;
        info.gmresInnerIterationsPerCycle = gmresInner;
        info.gmresOuterCycles = gmresOuter;
        info.operatorCalls = operatorCalls;
        info.preconditionerCalls = preconditionerCalls;
        info.updateCount = updateCount;
        info.completedIterations = lastCompleted;
        info.finalIterations = finalIterations;
        info.finalFlag = finalFlag;
        info.finalRelres = finalRelres;
        info.note = ['Live progress is inferred from Krylov operator and ', ...
            'preconditioner applications; the final completed count uses ', ...
            'RESVEC while finalIterations preserves MATLAB''s ITER output.'];
    end
end

function [maxWork, inner, outer, restarted] = ...
        krylovWorkLimits(solverName, maxIterations, restart, nTotal)
if strcmp(solverName, 'bicgstab')
    maxWork = maxIterations;
    inner = 1;
    outer = maxIterations;
    restarted = false;
    return
end

restarted = ~isempty(restart) && restart ~= nTotal;
if restarted
    inner = min(max(double(restart), 1), nTotal);
    outer = maxIterations;
else
    inner = min(maxIterations, nTotal);
    outer = 1;
end
maxWork = inner*outer;
end

function completed = gmresLiveIterations(preconditionerCalls, inner, restarted)
% The first preconditioner call belongs to the initial residual. Restarted
% GMRES adds one residual-preconditioning call after every inner cycle.
completed = max(preconditionerCalls - 1, 0);
if restarted
    completed = completed - floor(completed/(inner + 1));
end
end

function completed = exactCompletedIterations(solverName, resvec)
% RESVEC records every completed Krylov step even when ITER points to the
% best earlier approximation returned after a failed convergence attempt.
if strcmp(solverName, 'bicgstab')
    completed = 0.5*max(numel(resvec) - 1, 0);
else
    completed = max(numel(resvec) - 1, 0);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
function preconditioner = buildPreconditioner(sysInfo, params, solverName, ...
        nTotal, useGPU)
%BUILDPRECONDITIONER Build a matrix-free left preconditioner for Krylov.
%
% The basic Full-2D preconditioner is Jacobi:
%
%   M approx diag(A_diff + G_drift),  M^{-1}r = r ./ diag(M).
%
% This is cheap enough to be useful because the drift contribution is
% already diagonal and the DG/FV transport diagonal can be extracted from
% the Kronecker parts without assembling the global sparse matrix.
%
% As a stronger alternative, full2D_preconditioner='blockjacobi' builds one
% block per center-coordinate DG DOF. Each block contains all rho_x/rho_y
% unknowns for that center point, i.e. the contiguous vector segment
%
%   F(:, :, iX, iY).
%
% This drops center-neighbor DG couplings in M, but keeps the local sparse
% relative-coordinate transport, potential and CAP terms inside every block.

mode = preconditionerModeToString(readParam(params, ...
    'full2D_preconditioner', 'auto'));
maxDof = readParam(params, 'full2D_maxStoredPreconditionerDof', 2e7);

preconditioner = struct;
preconditioner.enabled = false;
preconditioner.apply = [];
preconditioner.info = struct;
preconditioner.info.requested = mode;
preconditioner.info.method = 'none';
preconditioner.info.enabled = false;
preconditioner.info.nTotal = nTotal;
preconditioner.info.maxStoredDof = maxDof;
preconditioner.info.executionDevice = ternary(useGPU, 'gpu', 'cpu');

if strcmp(solverName, 'direct')
    preconditioner.info.reason = 'Direct solver does not use a Krylov preconditioner.';
    return
end

switch mode
    case {'none', 'off', 'false', 'no'}
        preconditioner.info.reason = 'Preconditioner disabled by full2D_preconditioner.';
        return
    case {'rowabs', 'rowsum', 'row-sum', 'row-scaling'}
        preconditioner = enableRowAbsPreconditioner(preconditioner, ...
            sysInfo, params, nTotal, ...
            'Using row-absolute-sum scaling requested by full2D_preconditioner.');
        preconditioner = adaptPreconditionerForGPU(preconditioner, useGPU);
        return
    case {'blockjacobi', 'block-jacobi', 'rho-block-jacobi', ...
            'relative-block-jacobi', 'rho-blockjacobi'}
        if useGPU
            preconditioner = enableRowAbsPreconditioner(preconditioner, ...
                sysInfo, params, nTotal, ...
                ['GPU execution currently replaces the LU-based rho ', ...
                 'Block-Jacobi preconditioner by rowabs scaling.']);
            preconditioner.info.gpuFallbackFrom = 'blockjacobi-rho';
            preconditioner = adaptPreconditionerForGPU( ...
                preconditioner, true);
        else
            preconditioner = enableRelativeBlockJacobiPreconditioner( ...
                preconditioner, sysInfo, params, nTotal);
        end
        return
    case {'auto', 'jacobi', 'diagonal'}
        useJacobi = strcmp(mode, 'jacobi') || strcmp(mode, 'diagonal') ...
            || nTotal <= maxDof;
    otherwise
        error('DG:Full2D:UnknownPreconditioner', ...
            ['Unknown full2D_preconditioner "%s". Use auto, none, ', ...
             'jacobi, diagonal, rowabs or blockjacobi.'], mode);
end

if ~useJacobi
    preconditioner.info.reason = sprintf(['Jacobi preconditioner skipped for ', ...
        '%d DOFs because full2D_maxStoredPreconditionerDof=%d.'], ...
        nTotal, maxDof);
    return
end

allowLarge = readLogicalParam(params, 'full2D_allowLargePreconditioner', false);
if nTotal > maxDof && ~allowLarge
    preconditioner.info.reason = sprintf(['Jacobi preconditioner would store ', ...
        '%d diagonal entries, above full2D_maxStoredPreconditionerDof=%d. ', ...
        'Set full2D_allowLargePreconditioner=true only if that memory use is intended.'], ...
        nTotal, maxDof);
    return
end

if ~isfield(sysInfo, 'getDiagonal') || isempty(sysInfo.getDiagonal)
    preconditioner.info.reason = 'System diagonal is not available.';
    return
end

diagonal = sysInfo.getDiagonal();
[inverseDiagonal, floorInfo] = regularizePreconditionerDiagonal(diagonal, params);

fallbackOnSmallPivots = readLogicalParam(params, ...
    'full2D_jacobiFallbackOnSmallPivots', true);
if floorInfo.replacedEntries > 0 && fallbackOnSmallPivots
    preconditioner = enableRowAbsPreconditioner(preconditioner, ...
        sysInfo, params, nTotal, ...
        sprintf(['Requested Jacobi, but diag(A) contains %d small or zero ', ...
        'pivots; using row-absolute-sum scaling instead.'], ...
        floorInfo.replacedEntries));
    preconditioner.info.jacobiDiagonal = struct( ...
        'diagonalNnz', nnz(diagonal), ...
        'diagonalMinAbs', min(abs(diagonal)), ...
        'diagonalMaxAbs', max(abs(diagonal)), ...
        'floor', floorInfo);
    preconditioner = adaptPreconditionerForGPU(preconditioner, useGPU);
    return
end

preconditioner.enabled = true;
preconditioner.apply = @(r) inverseDiagonal.*r;
preconditioner.vector = inverseDiagonal;
preconditioner.info.method = 'jacobi';
preconditioner.info.enabled = true;
preconditioner.info.reason = ...
    'Using matrix-free Jacobi preconditioner with identity fallback for small pivots.';
preconditioner.info.diagonalNnz = nnz(diagonal);
preconditioner.info.diagonalMinAbs = min(abs(diagonal));
preconditioner.info.diagonalMaxAbs = max(abs(diagonal));
preconditioner.info.inverseDiagonalMinAbs = min(abs(inverseDiagonal));
preconditioner.info.inverseDiagonalMaxAbs = max(abs(inverseDiagonal));
preconditioner.info.floor = floorInfo;
preconditioner.info.note = ['Jacobi uses diag(A_diff)+diag(G_drift); it avoids ', ...
    'ILU-style global matrix assembly and is safe for the rho-basis ', ...
    'matrix-free path.'];
preconditioner = adaptPreconditionerForGPU(preconditioner, useGPU);
end

function tf = shouldRetryBlockJacobiWithRowAbs(flag, preconditioner, params)
tf = flag ~= 0 ...
    && preconditioner.enabled ...
    && isfield(preconditioner, 'info') ...
    && isfield(preconditioner.info, 'method') ...
    && strcmp(preconditioner.info.method, 'blockjacobi-rho') ...
    && readLogicalParam(params, ...
        'full2D_retryRowAbsOnBlockJacobiFailure', true);
end

function preconditioner = enableRowAbsPreconditioner(preconditioner, ...
        sysInfo, params, nTotal, reason)
%ENABLEROWABSPRECONDITIONER Use positive row-sum scaling as M^{-1}.
%
% The DG/FV Full-2D operator can have zero diagonal entries, so plain Jacobi
% may be singular or numerically harmful. Row-absolute-sum scaling is a
% sparse/matrix-free friendly fallback:
%
%   M_ii approx sum_j |A_ij|,    M^{-1}r = r ./ M_ii.
%
% It is less aggressive than diagonal Jacobi and does not require global
% matrix assembly.

if ~isfield(sysInfo, 'getRowAbsSum') || isempty(sysInfo.getRowAbsSum)
    preconditioner.info.reason = ...
        'Row-absolute-sum preconditioner requested, but row sums are not available.';
    return
end

rowAbsSum = sysInfo.getRowAbsSum();
[rowScale, scaleInfo] = regularizeRowAbsScale(rowAbsSum, params);

preconditioner.enabled = true;
preconditioner.apply = @(r) r./rowScale;
preconditioner.vector = rowScale;
preconditioner.info.method = 'rowabs';
preconditioner.info.enabled = true;
preconditioner.info.reason = reason;
preconditioner.info.nTotal = nTotal;
preconditioner.info.rowAbsMin = min(rowAbsSum);
preconditioner.info.rowAbsMax = max(rowAbsSum);
preconditioner.info.scaleMin = min(rowScale);
preconditioner.info.scaleMax = max(rowScale);
preconditioner.info.floor = scaleInfo;
preconditioner.info.note = ['Row-absolute-sum scaling approximates a ', ...
    'diagonal magnitude preconditioner without inverting zero diagonal ', ...
    'entries or assembling the full matrix.'];
end

function preconditioner = adaptPreconditionerForGPU(preconditioner, useGPU)
%ADAPTPRECONDITIONERFORGPU Transfer diagonal scaling data exactly once.
if ~useGPU || ~preconditioner.enabled
    return
end
if ~isfield(preconditioner, 'vector') || isempty(preconditioner.vector)
    error('DG:Full2D:UnsupportedGPUPreconditioner', ...
        'Preconditioner "%s" has no GPU scaling representation.', ...
        preconditioner.info.method);
end

gpuVector = gpuArray(preconditioner.vector);
switch preconditioner.info.method
    case 'jacobi'
        preconditioner.apply = @(r) gpuVector.*r;
    case 'rowabs'
        preconditioner.apply = @(r) r./gpuVector;
    otherwise
        error('DG:Full2D:UnsupportedGPUPreconditioner', ...
            ['Preconditioner "%s" is not supported on the GPU path. ', ...
             'Use rowabs or jacobi.'], preconditioner.info.method);
end
preconditioner = rmfield(preconditioner, 'vector');
preconditioner.info.executionDevice = 'gpu';
preconditioner.info.note = [preconditioner.info.note, ...
    ' Its scaling vector is stored on the selected GPU.'];
end

function preconditioner = enableRelativeBlockJacobiPreconditioner( ...
        preconditioner, sysInfo, params, nTotal)
%ENABLERELATIVEBLOCKJACOBIPRECONDITIONER Rho-block Jacobi for Full-2D.
%
% The global vector order is rho_x, rho_y, X-DG, Y-DG. Therefore all
% relative-coordinate DOFs for one center-coordinate DG point are contiguous
% in memory. The preconditioner assembles only the block diagonal
%
%   M = blockdiag(B_1, ..., B_nCenter),
%
% where B_c is the sparse rho_x/rho_y block of the full operator at fixed
% center DOF c. The full off-block DG couplings in X/Y are still applied by
% the Krylov operator A(u); they are only ignored in M^{-1}.

if ~isfield(sysInfo, 'getRelativeBlockData') ...
        || isempty(sysInfo.getRelativeBlockData)
    preconditioner.info.reason = ...
        'Block-Jacobi requested, but relative block data are not available.';
    return
end

nRelative = sysInfo.diff.relativeDof;
maxRelativeDof = readParam(params, ...
    'full2D_blockJacobiMaxRelativeDof', 512);
if nRelative > maxRelativeDof
    preconditioner.info.reason = sprintf(['Block-Jacobi skipped because ', ...
        'one rho block has %d DOFs, above full2D_blockJacobiMaxRelativeDof=%d.'], ...
        nRelative, maxRelativeDof);
    return
end

maxDof = readParam(params, 'full2D_maxBlockJacobiDof', 2e6);
allowLarge = readLogicalParam(params, 'full2D_allowLargePreconditioner', false);
if nTotal > maxDof && ~allowLarge
    preconditioner.info.reason = sprintf(['Block-Jacobi skipped because ', ...
        'the preconditioner would cover %d DOFs, above ', ...
        'full2D_maxBlockJacobiDof=%d. Set ', ...
        'full2D_allowLargePreconditioner=true only if this memory use is intended.'], ...
        nTotal, maxDof);
    return
end

storageMode = lower(char(readParam(params, ...
    'full2D_blockJacobiStorage', 'blockdiag')));
if ~any(strcmp(storageMode, {'blockdiag', 'explicit', 'sparse'}))
    error('DG:Full2D:UnknownBlockJacobiStorage', ...
        ['Unknown full2D_blockJacobiStorage "%s". The current ', ...
        'implementation supports "blockdiag".'], storageMode);
end

try
    blockData = sysInfo.getRelativeBlockData();
    [M, blockInfo] = assembleRelativeBlockJacobiMatrix(blockData, params);
    decompositionType = lower(char(readParam(params, ...
        'full2D_blockJacobiDecomposition', 'lu')));
    Mdec = decomposition(M, decompositionType);
catch err
    fallback = readLogicalParam(params, ...
        'full2D_blockJacobiFallbackToRowAbs', true);
    if fallback
        preconditioner = enableRowAbsPreconditioner(preconditioner, ...
            sysInfo, params, nTotal, ...
            sprintf(['Block-Jacobi setup failed (%s); using ', ...
            'row-absolute-sum scaling instead.'], err.message));
        preconditioner.info.blockJacobiError = err.message;
        return
    end
    rethrow(err)
end

preconditioner.enabled = true;
preconditioner.apply = @(r) Mdec\r;
preconditioner.info.method = 'blockjacobi-rho';
preconditioner.info.enabled = true;
preconditioner.info.reason = ['Using relative-coordinate Block-Jacobi ', ...
    'preconditioner requested by full2D_preconditioner.'];
preconditioner.info.nTotal = nTotal;
preconditioner.info.nBlocks = blockData.nCenter;
preconditioner.info.blockDof = blockData.nRelative;
preconditioner.info.blockMatrixNnz = nnz(M);
preconditioner.info.blockMatrixDensity = nnz(M)/numel(M);
preconditioner.info.storage = storageMode;
preconditioner.info.decomposition = decompositionType;
preconditioner.info.blockBuild = blockInfo;
preconditioner.info.note = ['Each preconditioner block acts on the local ', ...
    'rho_x/rho_y kernel for one fixed X/Y DG DOF. This is stronger than ', ...
    'scalar rowabs/Jacobi scaling but still avoids assembling the full ', ...
    'off-block DG transport matrix.'];
end

function [M, info] = assembleRelativeBlockJacobiMatrix(blockData, params)
%ASSEMBLERELATIVEBLOCKJACOBIMATRIX Explicit sparse block diagonal M.
%
% This matrix is not the full transport system. It stores only the
% center-block diagonal and can therefore be much sparser than A. The
% implementation intentionally starts with an explicit sparse blockdiag
% representation because it is easy to verify against small assembled
% systems and can later be replaced by a batched/matrix-free GPU apply.

nCenter = blockData.nCenter;
nRelative = blockData.nRelative;
blocks = cell(nCenter, 1);
blockNnz = zeros(nCenter, 1);
blockScale = zeros(nCenter, 1);
blockShift = zeros(nCenter, 1);
shift = readParam(params, 'full2D_blockJacobiShift', []);
relativeShift = readParam(params, 'full2D_blockJacobiRelativeShift', ...
    readParam(params, 'full2D_preconditionerFloor', 1e-10));
Irelative = speye(nRelative);

for ic = 1:nCenter
    block = sparse(blockData.getBlock(ic));
    if readLogicalParam(params, 'full2D_blockJacobiCheckSingular', true) ...
            && sprank(block) < nRelative
        error('DG:Full2D:SingularBlockJacobiBlock', ...
            ['Relative Block-Jacobi block %d has structural rank below ', ...
            'the block size %d. Use rowabs, choose an even rho grid, or ', ...
            'set full2D_blockJacobiCheckSingular=false for experiments.'], ...
            ic, nRelative);
    end
    blockScale(ic) = max(full(sum(abs(block), 2)));
    if isempty(shift)
        blockShift(ic) = relativeShift*max(blockScale(ic), 1);
    else
        blockShift(ic) = shift;
    end
    if blockShift(ic) ~= 0
        % The exact center-block diagonal can be singular because the
        % omitted neighboring DG couplings may regularize a local kernel.
        % A tiny shift is applied only to M, not to the physical system A.
        block = block + blockShift(ic)*Irelative;
    end
    blocks{ic} = block;
    blockNnz(ic) = nnz(block);
end

M = blkdiag(blocks{:});
M = sparse(M);

info = struct;
info.nBlocks = nCenter;
info.blockDof = nRelative;
info.totalDof = nCenter*nRelative;
info.shift = shift;
info.relativeShift = relativeShift;
info.shiftMode = ternary(isempty(shift), 'relative-to-local-rowabs', 'absolute');
info.blockScaleMin = min(blockScale);
info.blockScaleMax = max(blockScale);
info.blockShiftMin = min(blockShift);
info.blockShiftMax = max(blockShift);
info.blockNnzMin = min(blockNnz);
info.blockNnzMax = max(blockNnz);
info.blockNnzMean = mean(blockNnz);
info.note = ['B_c includes the local rho transport block and the local ', ...
    'sparse drift/CAP center block. Off-block X/Y DG couplings are omitted.'];
end

function value = ternary(condition, trueValue, falseValue)
if condition
    value = trueValue;
else
    value = falseValue;
end
end

function mode = preconditionerModeToString(rawMode)
if islogical(rawMode)
    if rawMode
        mode = 'jacobi';
    else
        mode = 'none';
    end
elseif isnumeric(rawMode)
    if rawMode ~= 0
        mode = 'jacobi';
    else
        mode = 'none';
    end
else
    mode = lower(char(rawMode));
end
end

function [inverseDiagonal, info] = regularizePreconditionerDiagonal(diagonal, params)
diagonal = diagonal(:);
maxAbs = max(abs(diagonal));
absoluteFloor = readParam(params, 'full2D_preconditionerAbsoluteFloor', []);
relativeFloor = readParam(params, 'full2D_preconditionerFloor', 1e-12);
if isempty(absoluteFloor)
    absoluteFloor = relativeFloor*max(maxAbs, 1);
end

% Do not replace zero pivots by a tiny floor and then divide by that value.
% For this DG/FV operator some diagonal entries are exactly zero; inverting a
% tiny artificial pivot amplifies those residual components and can trigger a
% BICGSTAB breakdown. Small-pivot entries are therefore left unpreconditioned
% by using the identity action on those components.
small = abs(diagonal) < absoluteFloor;
inverseDiagonal = complex(ones(size(diagonal)));
inverseDiagonal(~small) = 1./diagonal(~small);

info = struct;
info.absoluteFloor = absoluteFloor;
info.relativeFloor = relativeFloor;
info.replacedEntries = nnz(small);
info.smallPivotMode = 'identity';
info.note = ...
    'Small diagonal entries are not inverted; the preconditioner uses identity action there.';
end

function [rowScale, info] = regularizeRowAbsScale(rowAbsSum, params)
rowAbsSum = full(rowAbsSum(:));
maxAbs = max(rowAbsSum);
absoluteFloor = readParam(params, 'full2D_preconditionerRowAbsFloor', []);
if isempty(absoluteFloor)
    relativeFloor = readParam(params, 'full2D_preconditionerFloor', 1e-12);
    absoluteFloor = relativeFloor*max(maxAbs, 1);
else
    relativeFloor = NaN;
end

rowScale = rowAbsSum;
small = rowScale < absoluteFloor;
rowScale(small) = 1;

info = struct;
info.absoluteFloor = absoluteFloor;
info.relativeFloor = relativeFloor;
info.replacedEntries = nnz(small);
info.smallPivotMode = 'identity';
info.note = ...
    'Rows with near-zero absolute row sum are left unscaled.';
end

function value = readLogicalParam(params, name, defaultValue)
rawValue = readParam(params, name, defaultValue);
if islogical(rawValue)
    value = rawValue;
elseif isnumeric(rawValue)
    value = rawValue ~= 0;
else
    value = any(strcmpi(char(rawValue), {'true', 'yes', 'on', '1'}));
end
end

function [doSolve, status] = shouldSolve(mode, nTotal, maxAutoSolveDof)
if islogical(mode)
    doSolve = mode;
elseif isnumeric(mode)
    doSolve = mode ~= 0;
else
    switch lower(char(mode))
        case 'auto'
            doSolve = nTotal <= maxAutoSolveDof;
        case {'true', 'yes', 'on', 'solve', 'always'}
            doSolve = true;
        case {'false', 'no', 'off', 'operator-only', 'skip', 'never'}
            doSolve = false;
        otherwise
            error('DG:Full2D:UnknownSolveMode', ...
                'Unknown full2D_solve mode "%s".', char(mode));
    end
end

if doSolve
    status = 'solve-requested';
elseif ischar(mode) || isstring(mode)
    if strcmpi(char(mode), 'auto')
        status = 'operator-ready-solve-skipped-size';
    else
        status = 'operator-ready-solve-disabled';
    end
else
    status = 'operator-ready-solve-disabled';
end
end

function solverName = chooseSolver(params, matrixFreeOnly, nTotal, useGPU)
solverName = lower(char(readParam(params, 'full2D_solver', 'auto')));
if strcmp(solverName, 'auto')
    maxDirectSolveDof = readParam(params, 'full2D_maxDirectSolveDof', 150000);
    if useGPU
        solverName = 'bicgstab';
    elseif ~matrixFreeOnly && nTotal <= maxDirectSolveDof
        solverName = 'direct';
    else
        solverName = 'bicgstab';
    end
end
end

function value = gatherIfGPU(value)
if isa(value, 'gpuArray')
    value = gather(value);
end
end

function status = iterativeStatus(name, flag)
if flag == 0
    status = ['solved-', name];
else
    status = ['iterative-', name, '-not-converged'];
end
end

function tf = shouldStoreRHOArray(mat, nTotal)
params = getDGParams(mat);
mode = readParam(params, 'full2D_storeRHOArray', 'auto');
maxDof = readParam(params, 'full2D_maxStoredRHOArrayDof', 2e6);
if islogical(mode)
    tf = mode;
elseif isnumeric(mode)
    tf = mode ~= 0;
else
    switch lower(char(mode))
        case 'auto'
            tf = nTotal <= maxDof;
        case {'true', 'yes', 'on', 'store', 'always'}
            tf = true;
        case {'false', 'no', 'off', 'never'}
            tf = false;
        otherwise
            error('DG:Full2D:UnknownStorageMode', ...
                'Unknown full2D_storeRHOArray mode "%s".', char(mode));
    end
end
end

function params = getDGParams(mat)
params = struct;
if isfield(mat, 'dg') && isfield(mat.dg, 'params')
    params = mat.dg.params;
end
end

function value = readParam(params, name, defaultValue)
value = defaultValue;
if isstruct(params)
    if isfield(params, name) && ~isempty(params.(name))
        value = params.(name);
    elseif isfield(params, 'full2D') && isstruct(params.full2D) ...
            && isfield(params.full2D, name) && ~isempty(params.full2D.(name))
        value = params.full2D.(name);
    end
end
end

function text = modeToString(mode)
if islogical(mode)
    text = mat2str(mode);
elseif isnumeric(mode)
    text = num2str(mode);
else
    text = char(mode);
end
end
