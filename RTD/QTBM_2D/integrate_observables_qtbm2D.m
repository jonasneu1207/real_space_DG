function values = integrate_observables_qtbm2D(scattering,reservoir,~,~)
%INTEGRATE_OBSERVABLES_QTBM2D Spectral integrands per eV before quadrature.
% Rows: n (Nx), bond j (Nx-1), Landauer j, total incident flux magnitude.
% Positive current is positive electron flow from source to drain, as in DG.
h = scattering.hamiltonian; q = 1.602176634e-19; hbar = 1.054571817e-34;
[wL,wR] = contact_supply_qtbm2D(scattering.energyEV,reservoir);
L = scattering.sourceState; R = scattering.drainState;
n = (abs(L).^2.*wL+abs(R).^2.*wR)/(2*pi*h.dxM);
j = 2*q^2*h.hoppingEV/(2*pi*hbar)*( ...
    imag(conj(L(1:end-1,:)).*L(2:end,:)).*wL+ ...
    imag(conj(R(1:end-1,:)).*R(2:end,:)).*wR);
landauer = q^2/(2*pi*hbar)*scattering.transmission.*(wL-wR);
flux = q^2/(2*pi*hbar)*(wL+wR);
values = [n;j;landauer;flux];
end
