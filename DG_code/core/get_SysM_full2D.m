function [A, rhs, info] = get_SysM_full2D(mat, p, Vxy, EfL, EfR)
%GET_SYSM_FULL2D Boundary-focused system scaffold for full 2D Wigner-DG.
%
% No global 4D matrix is assembled in this first step. The function collects
% the separable components needed later:
%   get_Diff_full2D     rectangular DG/FV transport dimensions
%   get_Drift_full2D    relative-coordinate CAP scaffold
%   get_Boundary_full2D physical Source/Drain and specular Y boundaries

if nargin < 5
    EfR = p.EfR;
end
if nargin < 4
    EfL = p.EfL;
end

[A, rhs, diffInfo] = get_Diff_full2D(p);
[~, driftInfo] = get_Drift_full2D(mat, p, Vxy);
[boundary, boundaryInfo] = get_Boundary_full2D(mat, p, EfL, EfR, Vxy);

info = struct;
info.full2D = true;
info.fullMatrixAssembled = false;
info.diff = diffInfo;
info.drift = driftInfo;
info.boundary = boundaryInfo;
info.boundaryData = boundary;
info.dofOrder = p.index.order;
info.totalDofIfAssembled = p.index.nTotal;
end
