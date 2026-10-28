function poisson = prepare_poisson_qtbm2D(mat)
%PREPARE_POISSON_QTBM2D Assemble the shared RTD Poisson system without DG.
% X contacts: phi=0 and phi=bias. Y: periodic, final column is a duplicate.
rtdDir = fileparts(fileparts(mfilename('fullpath')));
oldPath = path;
restorePath = onCleanup(@() path(oldPath));
addpath(rtdDir);
get_RTD_doping_SI(mat);
[M,ids,types,info] = eval_Poisson2D_rtd(mat);
if ~strcmp(info.yBoundary,'periodic')
    error('QTBM2D:PoissonBoundary','The QTBM reference requires periodic Y.');
end
poisson = struct('matrix',M,'nodeIds',ids,'nodeTypes',types,'info',info,'mat',mat);
end
