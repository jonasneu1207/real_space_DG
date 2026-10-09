function [A, rhs, info] = get_SysM_full2D(mat, p, Vxy, EfL, EfR)
%GET_SYSM_FULL2D Full-2D rho-basis stationary system operator.
%
% The global unknown is
%
%   F = F(rho_x, rho_y, X, Y)
%
% with vector order
%
%   F(iRhoX, iRhoY, iX, iY) -> F(:),
%
% so rho_x is the fastest index, then rho_y, X-DG and Y-DG.
%
% Implemented system split:
%   Diff/transport: rectangular DG in X/Y with FV rho-operators and
%                   Source/Drain/specular boundary fluxes. The selectable
%                   position-dependent-mass path uses the full endpoint
%                   BenDaniel-Duke von-Neumann kinetic operator.
%   Drift/CAP:      FV-consistent sparse rho-basis potential blocks,
%                   projected locally over rectangular X/Y DG elements,
%                   plus a separable diagonal relative-coordinate CAP.
%
% For small systems this routine returns the assembled sparse matrix
%
%   A = A_diff + G_drift.
%
% For large systems A is left empty and info.apply(u) evaluates the same
% action matrix-free. The RHS is only materialized when the diff part is
% assembled, or when info.assembleRhs() is called explicitly.

if nargin < 5
    EfR = p.EfR;
end
if nargin < 4
    EfL = p.EfL;
end
if nargin < 3
    Vxy = [];
end

[boundary, boundaryInfo] = get_Boundary_full2D(mat, p, EfL, EfR, Vxy);
params = getDGParams(mat);
massModel = normalizeMassModel(readParam(params, ...
    'full2D_massModel', 'constant'));
switch massModel
    case 'constant'
        [A_diff, rhsDiff, diffInfo] = get_Diff_full2D(mat, p, boundary, false);
    case 'position-dependent-bdd'
        [A_diff, rhsDiff, diffInfo] = ...
            get_Diff_variableMass_full2D(mat, p, boundary, false);
end
[G_drift, driftInfo] = get_Drift_full2D(mat, p, Vxy);

fullMatrixAssembled = diffInfo.fullMatrixAssembled ...
    && driftInfo.fullPotentialAssembled;

if fullMatrixAssembled
    A = sparse(A_diff + G_drift);
    rhs = rhsDiff;
else
    A = [];
    rhs = [];
end

info = struct;
info.full2D = true;
info.massModel = massModel;
info.fullMatrixAssembled = fullMatrixAssembled;
info.matrixFreeAvailable = diffInfo.matrixFreeAvailable ...
    && driftInfo.matrixFreeAvailable;
info.reason = matrixReason(fullMatrixAssembled, diffInfo, driftInfo, p.index.nTotal);
info.diffMatrixAssembled = diffInfo.fullMatrixAssembled;
info.driftMatrixAssembled = driftInfo.fullPotentialAssembled;
info.diff = diffInfo;
info.drift = driftInfo;
info.boundary = boundaryInfo;
info.boundaryData = boundary;
info.dofOrder = p.index.order;
info.totalDofIfAssembled = p.index.nTotal;
info.apply = @(u) applySystem(diffInfo, driftInfo, u);
info.makeGPUApply = @() makeGPUSystemApply(diffInfo, driftInfo);
info.assemble = @() assembleSystem(diffInfo, driftInfo);
info.assembleRhs = @() diffInfo.assembleRhs();
info.getDiagonal = @() getSystemDiagonal(diffInfo, driftInfo);
info.getRowAbsSum = @() getSystemRowAbsSum(diffInfo, driftInfo);
info.getRelativeBlockData = @() getSystemRelativeBlockData(diffInfo, driftInfo);
if strcmp(get_YBoundaryType_full2D(mat), 'hard-wall')
    % Constrain the COMPLETE operator once. Constraining only transport
    % would let the potential/FV stencil regenerate forbidden entries.
    [A, rhs, info] = constrain_HardWall_full2D(A, rhs, info, p);
end
end

function model = normalizeMassModel(value)
model = lower(strrep(strtrim(char(value)), '_', '-'));
switch model
    case {'constant', 'scalar', 'legacy', 'reference'}
        model = 'constant';
    case {'position-dependent', 'spatial', 'variable', ...
            'variable-mass', 'bdd', 'ben-daniel-duke', ...
            'position-dependent-bdd'}
        model = 'position-dependent-bdd';
    otherwise
        error('DG:Full2D:UnknownMassModel', ...
            ['Unknown full2D_massModel "%s". Use "constant" or ', ...
             '"position-dependent".'], char(value));
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
if isstruct(params) && isfield(params, name) && ~isempty(params.(name))
    value = params.(name);
end
end

function y = applySystem(diffInfo, driftInfo, u)
%APPLYSYSTEM Matrix-free action of A_diff + G_drift.
y = diffInfo.apply(u) + driftInfo.apply(u);
end

function applyGPU = makeGPUSystemApply(diffInfo, driftInfo)
%MAKEGPUSYSTEMAPPLY Build one fully device-resident matrix-free operator.
if ~isfield(diffInfo, 'makeGPUApply') || isempty(diffInfo.makeGPUApply) ...
        || ~isfield(driftInfo, 'makeGPUApply') || isempty(driftInfo.makeGPUApply)
    error('DG:Full2D:MissingGPUOperator', ...
        'Diff and drift operators must both expose makeGPUApply().');
end
diffApplyGPU = diffInfo.makeGPUApply();
driftApplyGPU = driftInfo.makeGPUApply();
applyGPU = @(u) diffApplyGPU(u) + driftApplyGPU(u);
end

function A = assembleSystem(diffInfo, driftInfo)
%ASSEMBLESYSTEM Explicit sparse assembly for small diagnostic problems.
A = sparse(diffInfo.assemble() + driftInfo.assemble());
end

function diagonal = getSystemDiagonal(diffInfo, driftInfo)
%GETSYSTEMDIAGONAL Diagonal of A_diff + G_drift for Jacobi preconditioning.
diagonal = diffInfo.getDiagonal() + driftInfo.getDiagonal();
end

function rowAbsSum = getSystemRowAbsSum(diffInfo, driftInfo)
%GETSYSTEMROWABSSUM Row magnitude scaling for matrix-free preconditioning.
rowAbsSum = diffInfo.getRowAbsSum() + driftInfo.getRowAbsSum();
end

function blockData = getSystemRelativeBlockData(diffInfo, driftInfo)
%GETSYSTEMRELATIVEBLOCKDATA Block-Jacobi data for rho_x/rho_y kernels.
%
% The Full-2D vector is ordered as
%
%   F(iRhoX,iRhoY,iX,iY) -> F(:),
%
% so the contiguous entries
%
%   ((c-1)*nRelative+1) : c*nRelative
%
% contain all relative-coordinate unknowns for one fixed center-coordinate
% DG DOF c. The Block-Jacobi preconditioner uses exactly these contiguous
% rho-blocks:
%
%   B_c = A_diff(c,c in center space) + G_drift(c,c in center space).
%
% This keeps the preconditioner local in X/Y, but it retains the sparse
% rho_x/rho_y transport coupling and the local sparse potential/CAP term.

if ~isfield(diffInfo, 'getRelativeBlockData') ...
        || isempty(diffInfo.getRelativeBlockData)
    error('DG:Full2D:MissingBlockData', ...
        'Diff operator does not expose relative block data.');
end

diffBlockData = diffInfo.getRelativeBlockData();
if isfield(driftInfo, 'getRelativeBlockData') ...
        && ~isempty(driftInfo.getRelativeBlockData)
    driftBlockData = driftInfo.getRelativeBlockData();
else
    driftDiagonal = reshape(driftInfo.getDiagonal(), ...
        diffBlockData.nRelative, diffBlockData.nCenter);
    driftBlockData = struct;
    driftBlockData.getBlock = @(centerId) spdiags( ...
        driftDiagonal(:, centerId), 0, diffBlockData.nRelative, ...
        diffBlockData.nRelative);
end

blockData = struct;
blockData.nCenter = diffBlockData.nCenter;
blockData.nRelative = diffBlockData.nRelative;
blockData.nTotal = diffBlockData.nTotal;
blockData.diff = diffBlockData;
blockData.drift = driftBlockData;
blockData.blockType = 'relative-rho-blocks-per-center-dof';
blockData.dofOrder = {'rho_x', 'rho_y', 'X-DG-DOF', 'Y-DG-DOF'};
blockData.note = ['System block c is the local rho_x/rho_y transport ', ...
    'block plus the sparse drift/CAP center block for center DOF c. ', ...
    'Off-block center-coordinate DG couplings remain in A, not in M.'];
blockData.getBlock = @(centerId) systemRelativeBlock(diffBlockData, ...
    driftBlockData, centerId);
blockData.getCoupledBlock = @(rowId, colId) ...
    diffBlockData.getCoupledBlock(rowId,colId) ...
    + driftBlockData.getCoupledBlock(rowId,colId);
end

function block = systemRelativeBlock(diffBlockData, driftBlockData, centerId)
block = diffBlockData.getBlock(centerId) ...
    + driftBlockData.getBlock(centerId);
block = sparse(block);
end

function reason = matrixReason(assembled, diffInfo, driftInfo, nTotal)
if assembled
    reason = sprintf('Full-2D system matrix assembled for %d total DOFs.', nTotal);
else
    reason = sprintf(['Full-2D system kept matrix-free for %d total DOFs. ', ...
        'Diff: %s Drift: %s'], nTotal, diffInfo.reason, driftInfo.reason);
end
end
