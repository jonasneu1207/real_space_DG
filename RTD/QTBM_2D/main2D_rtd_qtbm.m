function result = main2D_rtd_qtbm(input)
%MAIN2D_RTD_QTBM Independent Y-homogeneous periodic QTBM RTD simulation.
% result = main2D_rtd_qtbm(struct('bias',0)); % self-consistent equilibrium
% result = main2D_rtd_qtbm(struct('bias',0:0.02:0.6)); % IV sweep
% result = main2D_rtd_qtbm(struct('bias',0,'poissonXBoundary','neumann'));
% No clearing of workspace, plots or automatic writes.
if nargin<1, input = struct; end
result = run_RTD_qtbm2D(input);
end
