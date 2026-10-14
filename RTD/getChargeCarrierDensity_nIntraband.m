function n = getChargeCarrierDensity_nIntraband(E, V, muf, beta, m)
% function getChargeCarrierDensity(E, V, muf, beta, m) evaluates the
% electron density for bulk semiconductor

h  = 6.626E-34;
hb = h/2/pi;
q  = 1.602E-19;
m0 = 9.109E-31;

dE = E(2)-E(1);
% density of states 
Nc  = (2*m*m0)^(3/2)/2/pi^2/hb^3;
Dc  = Nc*sqrt(q*(E-V));
% disp("q*(E-V): ");
% disp(sqrt(q*(E-V)));
n = zeros(1, length(muf));

% evaluate the charge carrier for each fermi level
for IF = 1 : length(muf)
    n(IF) = sum(Dc./(1+exp(beta*(E-muf(IF)))))*dE*q;
end