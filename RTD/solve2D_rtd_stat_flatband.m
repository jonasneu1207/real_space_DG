function results = solve2D_rtd_stat_flatband(solver,mat)
%SOLVE2D_RTD_STAT_FLATBAND Fixed-potential RTD sweep using the shared solver.
% No Poisson iteration: a prescribed voltage ramp is added between leads.
if ~strcmpi(solver,'DG')
    error('RTD:UnsupportedSolver','This RTD setup uses solve_transport_DG_full2D.');
end
bias = mat.Vdd(:).';
x = mat.x(:);
rightLeadStart = x(end)-mat.Lc2;
ramp = min(1,max(0,(x-mat.Lc1)/(rightLeadStart-mat.Lc1)));
sweep = struct('bias',bias);
for ib = 1:numel(bias)
    V = mat.Vc-bias(ib)*ramp;
    efL = mat.Vc(1,1)+0.05;
    efR = mat.Vc(end,1)+0.05-bias(ib);
    if isfield(mat,'EfL'), efL = mat.EfL; end
    if isfield(mat,'EfR'), efR = mat.EfR-bias(ib); end
    results = solve_transport_DG_full2D(mat,V,efL,efR);
    if isempty(results.rho)
        error('RTD:SolveSkipped','Enable full2D_solve to compute an RTD state.');
    end
    sweep.current(:,ib) = results.j;
    sweep.density(:,:,ib) = results.n;
    sweep.potential(:,:,ib) = V;
    sweep.flag(ib) = results.info.solve.flag;
    sweep.relres(ib) = results.info.solve.relres;
    sweep.periodicYStateError(ib) = results.info.solve.periodicYStateError;
    if sweep.flag(ib) ~= 0
        warning('RTD:NotConverged', ...
            'Bias %g V: flag=%g, true relative residual=%g. Do not interpret as converged.', ...
            bias(ib),sweep.flag(ib),sweep.relres(ib));
    end
end
% Keep the full state only for the last bias; retain observables for the sweep.
results.sweep = sweep;
results.Id = mean(sweep.current,1);
results.n_electron = sweep.density;
results.Vc = sweep.potential;
end
