function [integrated,info,spectrum] = adaptive_energy_qtbm2D(h,reservoir,mat,options)
%ADAPTIVE_ENERGY_QTBM2D Adaptive 8/16-point quadrature with band-edge splits.
% A sin^2 map removes the integrable inverse-square-root edge singularity.
% energyStepEV bounds initial panel width; it must also resolve resonances.
edges = [h.sourceEdgeEV,h.drainEdgeEV];
lower = min(edges);
upper = max([reservoir.source.muEV-min(reservoir.source.transverseEnergyEV), ...
    reservoir.drain.muEV-min(reservoir.drain.transverseEnergyEV),max(edges)]) ...
    +options.energyTailKT*reservoir.kTEV;
if upper>=min(edges)+4*h.hoppingEV
    error('QTBM2D:BandLimit','Refine dx: the occupied window reaches the artificial FD band edge.');
end
cuts = unique([lower,edges(edges>lower & edges<upper),upper]);
a = []; b = [];
for k = 1:numel(cuts)-1
    grid = linspace(cuts(k),cuts(k+1),max(1,ceil((cuts(k+1)-cuts(k))/options.energyStepEV))+1);
    a = [a,grid(1:end-1)]; b = [b,grid(2:end)]; %#ok<AGROW>
end
N = h.size; integrated = zeros(2*N+1,1); errors = integrated;
if 24*numel(a)>options.maxEnergyEvaluations
    error('QTBM2D:EnergyBudget','Initial energy grid exceeds maxEnergyEvaluations.');
end
evaluations = 0; totalPanels = 0; spectrumE = []; spectrumT = [];
scale = []; converged = false;
for level = 0:options.maxEnergyRefinements
    high = zeros(2*N+1,numel(a)); low = high;
    allE = zeros(16,numel(a)); allT = allE;
    for first = 1:16:numel(a)
        ids = first:min(first+15,numel(a));
        for order = [8,16]
            [nodes,weights] = gauss_rule_qtbm2D(order);
            theta = (nodes+1)*pi/4;
            E = a(ids)+(b(ids)-a(ids)).*sin(theta).^2;
            W = weights*pi/4.*(b(ids)-a(ids)).*sin(2*theta);
            states = solve_scattering_qtbm2D(h,E(:));
            values = integrate_observables_qtbm2D(states,reservoir,mat,options);
            values = reshape(values,2*N+1,order,numel(ids));
            moments = reshape(sum(values.*reshape(W,1,order,numel(ids)),2),2*N+1,numel(ids));
            if order==16
                high(:,ids) = moments; allE(:,ids) = E;
                allT(:,ids) = reshape(states.transmission,16,numel(ids));
            else
                low(:,ids) = moments;
            end
            evaluations = evaluations+numel(E);
        end
    end
    if any(~isfinite(high(:))) || any(~isfinite(low(:)))
        error('QTBM2D:ScatteringFailure','Nonfinite spectral integrals; inspect the potential and grid.');
    end
    delta = abs(high-low);
    if isempty(scale)
        estimate = sum(high,2);
        nScale = max([max(abs(estimate(1:N))),max(mat.Nd(:)),1]);
        jScale = max([abs(estimate(2*N)),estimate(end)*1e-8,1]);
        scale = [repmat(nScale,N,1);repmat(jScale,N,1);max(estimate(end),1)];
    end
    budget = options.energyRelativeTolerance*scale.*((b-a)/(upper-lower));
    accepted = all(delta<=budget,1);
    if level==options.maxEnergyRefinements, accepted(:) = true; end
    integrated = integrated+sum(high(:,accepted),2);
    errors = errors+sum(delta(:,accepted),2);
    spectrumE = [spectrumE;reshape(allE(:,accepted),[],1)]; %#ok<AGROW>
    spectrumT = [spectrumT;reshape(allT(:,accepted),[],1)]; %#ok<AGROW>
    totalPanels = totalPanels+sum(accepted);
    if all(accepted)
        converged = all(errors<=options.energyRelativeTolerance*scale);
        break
    end
    a = a(~accepted); b = b(~accepted); mid = (a+b)/2;
    a = reshape([a;mid],1,[]); b = reshape([mid;b],1,[]);
    if evaluations+24*numel(a)>options.maxEnergyEvaluations
        error('QTBM2D:EnergyBudget','Energy quadrature exceeded maxEnergyEvaluations; no result accepted.');
    end
end
[spectrumE,order] = sort(spectrumE);
spectrum = struct('energyEV',spectrumE,'transmission',spectrumT(order));
info = struct('converged',converged,'evaluations',evaluations,'panels',totalPanels, ...
    'refinementLevels',level,'energyWindowEV',[lower,upper], ...
    'densityErrorEstimate',max(errors(1:N))/scale(1), ...
    'currentErrorEstimate',max(errors(N+1:2*N))/scale(N+1), ...
    'tailKT',options.energyTailKT,'relativeTolerance',options.energyRelativeTolerance);
if ~converged
    error('QTBM2D:EnergyConvergence','Energy quadrature failed: density %.3g, current %.3g.', ...
        info.densityErrorEstimate,info.currentErrorEstimate);
end
end
