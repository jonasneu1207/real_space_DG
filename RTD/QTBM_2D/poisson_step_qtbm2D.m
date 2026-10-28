function [phiCandidate,info] = poisson_step_qtbm2D(state,poisson)
%POISSON_STEP_QTBM2D Shared direct/Newton step, independent of DG data types.
% state.phi: [Nx,Ny] electrostatic potential [V]. state.bias: drain bias [V].
% state.nPoisson: [Nx,Ny] CONTROL-VOLUME average density [m^-3], including
% the duplicate final Y column. Point samples alone do not meet this contract.
% info.converged only describes the linear update, NOT self-consistency.
rtdDir = fileparts(fileparts(mfilename('fullpath')));
oldPath = path;
restorePath = onCleanup(@() path(oldPath));
addpath(rtdDir);
[phiCandidate,~,info] = solve_Poisson2D_rtd(poisson.matrix, ...
    poisson.nodeTypes,poisson.nodeIds,state.nPoisson,poisson.mat,state.phi,state.bias);
end
