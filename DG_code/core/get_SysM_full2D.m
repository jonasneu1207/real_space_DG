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
%                   Source/Drain/specular boundary fluxes.
%   Drift/CAP:      diagonal rho-basis potential difference plus separable
%                   relative-coordinate CAP.
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
[A_diff, rhsDiff, diffInfo] = get_Diff_full2D(mat, p, boundary);
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
info.assemble = @() assembleSystem(diffInfo, driftInfo);
info.assembleRhs = @() diffInfo.assembleRhs();
info.getDiagonal = @() getSystemDiagonal(diffInfo, driftInfo);
info.getRowAbsSum = @() getSystemRowAbsSum(diffInfo, driftInfo);
end

function y = applySystem(diffInfo, driftInfo, u)
%APPLYSYSTEM Matrix-free action of A_diff + G_drift.
y = diffInfo.apply(u) + driftInfo.apply(u);
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

function reason = matrixReason(assembled, diffInfo, driftInfo, nTotal)
if assembled
    reason = sprintf('Full-2D system matrix assembled for %d total DOFs.', nTotal);
else
    reason = sprintf(['Full-2D system kept matrix-free for %d total DOFs. ', ...
        'Diff: %s Drift: %s'], nTotal, diffInfo.reason, driftInfo.reason);
end
end
