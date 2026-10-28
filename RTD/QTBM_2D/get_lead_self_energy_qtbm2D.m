function lead = get_lead_self_energy_qtbm2D(energyEV,hamiltonian,side)
%GET_LEAD_SELF_ENERGY_QTBM2D Exact retarded lead for the discrete Hx.
% injection=sqrt(Gamma) gives a spectral state G*sqrt(Gamma) [eV^-1/2].
side = validatestring(side,{'source','drain'});
t = hamiltonian.hoppingEV; edge = hamiltonian.([side,'EdgeEV']);
kinetic = energyEV-edge;
inside = kinetic>0 & kinetic<4*t;
sigma = complex(zeros(size(kinetic)));
% Product form preserves accuracy close to a propagation threshold.
sigma(inside) = kinetic(inside)/2-t ...
    -0.5i*sqrt(kinetic(inside).* (4*t-kinetic(inside)));
outside = ~inside;
% Reciprocal form avoids cancellation far outside the FD band.
u = kinetic(outside)-2*t;
root = sqrt(kinetic(outside).*(kinetic(outside)-4*t));
sigma(outside) = 2*t^2./(u+sign(u).*root);
gamma = max(0,-2*imag(sigma));
lead = struct('selfEnergyEV',sigma,'gammaEV',gamma, ...
    'injection',sqrt(gamma),'propagating',inside, ...
    'velocityMS',gamma*1.602176634e-19*hamiltonian.dxM/1.054571817e-34);
end
