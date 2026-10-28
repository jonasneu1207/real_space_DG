function scattering = solve_scattering_qtbm2D(hamiltonian,energyEV)
%SOLVE_SCATTERING_QTBM2D Two retarded spectral states per energy.
% Batched tridiagonal elimination; each column is an independent energy.
E = energyEV(:).'; N = hamiltonian.size; t = hamiltonian.hoppingEV;
left = get_lead_self_energy_qtbm2D(E,hamiltonian,'source');
right = get_lead_self_energy_qtbm2D(E,hamiltonian,'drain');
d = complex(E-hamiltonian.diagonalEV);
d(1,:) = d(1,:)-left.selfEnergyEV;
d(end,:) = d(end,:)-right.selfEnergyEV;
psiL = complex(zeros(N,numel(E))); psiR = psiL;
psiL(1,:) = left.injection; psiR(end,:) = right.injection;
for k = 2:N
    factor = t./d(k-1,:);
    d(k,:) = d(k,:)-t*factor;
    psiL(k,:) = psiL(k,:)-factor.*psiL(k-1,:);
end
psiL(end,:) = psiL(end,:)./d(end,:);
psiR(end,:) = psiR(end,:)./d(end,:);
for k = N-1:-1:1
    psiL(k,:) = (psiL(k,:)-t*psiL(k+1,:))./d(k,:);
    psiR(k,:) = -t*psiR(k+1,:)./d(k,:);
end
scattering = struct('energyEV',E,'sourceState',psiL,'drainState',psiR, ...
    'transmission',right.gammaEV.*abs(psiL(end,:)).^2, ...
    'reverseTransmission',left.gammaEV.*abs(psiR(1,:)).^2, ...
    'hamiltonian',hamiltonian);
end
