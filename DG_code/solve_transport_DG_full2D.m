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
%   - local diagonal potential difference in rho basis,
%   - separable CAP on the outer relative-coordinate boundaries.
%
% Matrix policy:
%   Small systems are assembled as sparse matrices and solved directly by
%   default. Large systems keep A empty and expose sysInfo.apply(u). A
%   matrix-free Krylov solve is available via params.full2D_solve = true
%   and params.full2D_solver = 'bicgstab' or 'gmres'.

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
maxAutoSolveDof = readParam(params, 'full2D_maxAutoSolveDof', 6e6);
maxMatrixFreeSolveDof = readParam(params, 'full2D_maxMatrixFreeSolveDof', 2e7);

rho = [];
solveInfo = struct;
solveInfo.requestedMode = mode;
solveInfo.nTotal = nTotal;
solveInfo.matrixAssembled = sysInfo.fullMatrixAssembled;
solveInfo.matrixFreeAvailable = sysInfo.matrixFreeAvailable;
solveInfo.tolerance = readParam(params, 'full2D_solverTol', 1e-6);
solveInfo.maxIterations = readParam(params, 'full2D_solverMaxIt', 200);

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

solverName = chooseSolver(params, isempty(A), nTotal);
solveInfo.solver = solverName;
preconditioner = buildPreconditioner(sysInfo, params, solverName, nTotal);
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
        applyA = @(u) sysInfo.apply(u);
        if preconditioner.enabled
            [rho, flag, relres, iter, resvec] = bicgstab(applyA, rhs, ...
                solveInfo.tolerance, solveInfo.maxIterations, ...
                preconditioner.apply);
        else
            [rho, flag, relres, iter, resvec] = bicgstab(applyA, rhs, ...
                solveInfo.tolerance, solveInfo.maxIterations);
        end
        solveInfo.status = iterativeStatus('bicgstab', flag);
        solveInfo.flag = flag;
        solveInfo.relres = relres;
        solveInfo.iterations = iter;
        solveInfo.resvec = resvec;

    case 'gmres'
        applyA = @(u) sysInfo.apply(u);
        restart = readParam(params, 'full2D_gmresRestart', []);
        if preconditioner.enabled
            [rho, flag, relres, iter, resvec] = gmres(applyA, rhs, restart, ...
                solveInfo.tolerance, solveInfo.maxIterations, ...
                preconditioner.apply);
        else
            [rho, flag, relres, iter, resvec] = gmres(applyA, rhs, restart, ...
                solveInfo.tolerance, solveInfo.maxIterations);
        end
        solveInfo.status = iterativeStatus('gmres', flag);
        solveInfo.flag = flag;
        solveInfo.relres = relres;
        solveInfo.iterations = iter;
        solveInfo.resvec = resvec;

    otherwise
        error('DG:Full2D:UnknownSolver', ...
            'Unknown full2D_solver "%s". Use auto, direct, bicgstab or gmres.', ...
            solverName);
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
function preconditioner = buildPreconditioner(sysInfo, params, solverName, nTotal)
%BUILDPRECONDITIONER Build a matrix-free left preconditioner for Krylov.
%
% The first Full-2D preconditioner is Jacobi:
%
%   M approx diag(A_diff + G_drift),  M^{-1}r = r ./ diag(M).
%
% This is cheap enough to be useful because the drift contribution is
% already diagonal and the DG/FV transport diagonal can be extracted from
% the Kronecker parts without assembling the global sparse matrix.

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
        return
    case {'auto', 'jacobi', 'diagonal'}
        useJacobi = strcmp(mode, 'jacobi') || strcmp(mode, 'diagonal') ...
            || nTotal <= maxDof;
    otherwise
        error('DG:Full2D:UnknownPreconditioner', ...
            ['Unknown full2D_preconditioner "%s". Use auto, none, ', ...
             'jacobi, diagonal or rowabs.'], mode);
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
    return
end

preconditioner.enabled = true;
preconditioner.apply = @(r) inverseDiagonal.*r;
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

function solverName = chooseSolver(params, matrixFreeOnly, nTotal)
solverName = lower(char(readParam(params, 'full2D_solver', 'auto')));
if strcmp(solverName, 'auto')
    maxDirectSolveDof = readParam(params, 'full2D_maxDirectSolveDof', 150000);
    if ~matrixFreeOnly && nTotal <= maxDirectSolveDof
        solverName = 'direct';
    else
        solverName = 'bicgstab';
    end
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
