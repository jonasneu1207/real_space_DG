function result = run_RTD_qtbm2D(input)
%RUN_RTD_QTBM2D Prepare, fixed-potential sweep, or independent Poisson loop.
% Returns result.points, result.sweep, options and mat. No automatic saves.
if nargin<1, input = struct; end
options = initParams_qtbm2D(input);
rtdDir = fileparts(fileparts(mfilename('fullpath')));
oldPath = path; restorePath = onCleanup(@() path(oldPath));
addpath(rtdDir);
mat = buildDevice2D_rtd('n-RTD-constM-comp',options.Ly,options.dy);
mat.Vdd = options.bias;
mat.poisson.yBoundary = 'periodic';
mat.poisson.opt.solve = options.poissonMethod; mat.poisson.opt.alpha = 1;
ramp = bias_ramp_qtbm2D(mat);
phi0 = repmat(options.bias(1)*ramp,1,mat.Ny);
potentialEV = mat.Vc-phi0;
contacts = get_contacts_qtbm2D(mat,potentialEV,options);
poisson = prepare_poisson_qtbm2D(mat);
if strcmp(options.simulation,'prepare')
    result = struct('status','prepared','transportImplemented',true, ...
        'options',options,'mat',mat,'phi0',phi0,'potentialEV',potentialEV, ...
        'contacts',contacts,'poisson',poisson);
    return
end
if strcmp(options.simulation,'self-consistent')
    result = solve_self_consistent_qtbm2D(mat,poisson,options);
else
    points = cell(1,numel(options.bias));
    for k = 1:numel(options.bias)
        phi = repmat(options.bias(k)*ramp,1,mat.Ny); potentialEV = mat.Vc-phi;
        contacts = get_contacts_qtbm2D(mat,potentialEV,options);
        transport = solve_transport_qtbm2D(mat,potentialEV,contacts,options);
        [nPoisson,projection] = project_density_qtbm2D(transport,mat);
        points{k} = struct('bias',options.bias(k),'phi',phi,'potentialEV',potentialEV, ...
            'transport',transport,'nPoisson',nPoisson,'projection',projection, ...
            'converged',true,'iterations',0);
        if options.verbose
            fprintf('QTBM fixed V=%g V: j=%g A/m^2\n',options.bias(k),transport.currentDensity);
        end
    end
    result = collect_sweep_qtbm2D(points,mat,options,'fixed-potential-complete');
end
end
