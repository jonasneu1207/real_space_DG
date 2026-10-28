function bound = get_bound_states_qtbm2D(h,reservoir,options)
%GET_BOUND_STATES_QTBM2D True states below BOTH lead continua.
% Evanescent tails enter normalization. Equilibrium occupation is unique;
% biased isolated states require an explicit population assumption.
upper = min(h.sourceEdgeEV,h.drainEdgeEV)-1e-10;
lower = min(h.potentialEV)-1;
count = countBelow(upper,h);
bound = struct('energyEV',zeros(0,1),'states',zeros(h.size,0), ...
    'density',zeros(h.size,1),'occupation','none');
if count==0, return; end
equilibrium = abs(reservoir.source.muEV-reservoir.drain.muEV)<1e-9;
if equilibrium
    population = 'equilibrium';
else
    population = options.boundStateOccupation;
    if strcmp(population,'error')
        error('QTBM2D:BoundStateOccupation', ...
            ['True bound states below both contacts have no uniquely determined biased occupation. ', ...
            'Set boundStateOccupation explicitly to source, drain, or average, and check its influence.']);
    end
end
bound.occupation = population;
for k = 1:count
    lo = lower; hi = upper;
    for iteration = 1:80
        mid = (lo+hi)/2;
        if countBelow(mid,h)>=k, hi = mid; else, lo = mid; end
        if hi-lo<1e-12, break; end
    end
    E = (lo+hi)/2;
    left = get_lead_self_energy_qtbm2D(E,h,'source');
    right = get_lead_self_energy_qtbm2D(E,h,'drain');
    H = h.matrixEV;
    H(1,1) = H(1,1)+left.selfEnergyEV;
    H(end,end) = H(end,end)+right.selfEnergyEV;
    [psi,~,flag] = eigs(real(H),1,E+1e-7,struct('tol',1e-11));
    if flag~=0, error('QTBM2D:BoundEigenvector','Bound-state eigenvector did not converge.'); end
    lambdaL = abs(left.selfEnergyEV/h.hoppingEV);
    lambdaR = abs(right.selfEnergyEV/h.hoppingEV);
    norm2 = sum(abs(psi).^2)+abs(psi(1))^2*lambdaL^2/(1-lambdaL^2) ...
        +abs(psi(end))^2*lambdaR^2/(1-lambdaR^2);
    psi = psi/sqrt(norm2);
    [wL,wR] = contact_supply_qtbm2D(E,reservoir);
    switch population
        case {'equilibrium','source'}, weight = wL;
        case 'drain', weight = wR;
        case 'average', weight = (wL+wR)/2;
    end
    bound.energyEV(end+1,1) = E;
    bound.states(:,end+1) = psi;
    bound.density = bound.density+abs(psi).^2/h.dxM*weight;
end
end

function count = countBelow(E,h)
% Inertia of H+Sigma(E)-E; eigenvalues decrease monotonically below bands.
left = get_lead_self_energy_qtbm2D(E,h,'source');
right = get_lead_self_energy_qtbm2D(E,h,'drain');
d = h.diagonalEV-E;
d(1) = d(1)+real(left.selfEnergyEV); d(end) = d(end)+real(right.selfEnergyEV);
pivot = d(1); count = double(pivot<0);
for k = 2:h.size
    if abs(pivot)<1e-100, pivot = -1e-100; end
    pivot = d(k)-h.hoppingEV^2/pivot;
    count = count+double(pivot<0);
end
end
