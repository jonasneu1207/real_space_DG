function transport = solve_transport_qtbm2D(mat,potentialEV,contacts,options)
%SOLVE_TRANSPORT_QTBM2D Separable periodic-Y open quantum transport.
% n: point density [m^-3], jx: bond current density [A/m^2].
% currentPerWidth and j: Y-integrated current [A/m], as in the DG solver.
h = build_hamiltonian_qtbm2D(mat,potentialEV);
reservoir = get_occupation_qtbm2D(mat,contacts,[],options);
probe = linspace(min(potentialEV(:)),max(reservoir.source.muEV,reservoir.drain.muEV)+ ...
    10*reservoir.kTEV,31);
[wL,wR] = contact_supply_qtbm2D(probe,reservoir);
fine = reservoir; fine.kzPoints = 2*reservoir.kzPoints;
[vL,vR] = contact_supply_qtbm2D(probe,fine);
kzError = max(abs([vL-wL,vR-wR]))/max([vL,vR]);
if kzError>options.energyRelativeTolerance/4
    error('QTBM2D:KzConvergence','Increase kzPoints: kz occupation error estimate %.3g.',kzError);
end
[integrated,energyInfo,spectrum] = adaptive_energy_qtbm2D(h,reservoir,mat,options);
if min(spectrum.transmission)<0 || max(spectrum.transmission)>1+options.currentTolerance
    error('QTBM2D:ScatteringUnitarity','Transmission lies outside the physical interval [0,1].');
end
bound = get_bound_states_qtbm2D(h,reservoir,options);
nx = h.size; density = integrated(1:nx)+bound.density;
jx = integrated(nx+1:2*nx-1); landauer = integrated(2*nx);
currentScale = max([abs(landauer),integrated(end)*1e-8,1]);
currentError = max(abs(jx-landauer))/currentScale;
if currentError>options.currentTolerance
    error('QTBM2D:CurrentConservation','Bond and Landauer currents disagree: %.3g.',currentError);
end
transport = struct('status','converged','converged',true,'potentialEV',potentialEV, ...
    'n',repmat(density,1,mat.Ny),'jx',repmat(jx,1,mat.Ny), ...
    'currentDensity',landauer,'currentPerWidth',jx*contacts.widthM, ...
    'j',landauer*contacts.widthM,'xBondM',(h.xM(1:end-1)+h.xM(2:end))/2, ...
    'reservoir',reservoir,'contacts',contacts,'spectrum',spectrum,'boundStates',bound, ...
    'diagnostics',struct('energy',energyInfo,'kzRelativeError',kzError, ...
    'currentConservationError',currentError,'densityMinimum',min(density)));
end
