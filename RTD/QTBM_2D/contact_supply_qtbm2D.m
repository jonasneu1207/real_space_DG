function [source,drain] = contact_supply_qtbm2D(energyEV,reservoir)
%CONTACT_SUPPLY_QTBM2D Sum discrete ky; integrate kz/(2*pi) [m^-2].
% Includes spin 2 and deg_factor, divided by ONE periodic Y width.
E = energyEV(:).'; kT = reservoir.kTEV;
[z,w] = gauss_rule_qtbm2D(reservoir.kzPoints);
values = cell(1,2); sides = {'source','drain'};
for s = 1:2
    c = reservoir.(sides{s});
    eta = (c.muEV-E-c.transverseEnergyEV(:))/kT;
    umax = sqrt(max(40,max(eta(:))+40));
    u = (z+1)*umax/2; wu = w*umax/2;
    value = zeros(size(E));
    for mode = 1:size(eta,1)
        f = 1./(1+exp(min(700,u.^2-eta(mode,:))));
        value = value+wu.'*f;
    end
    kzScale = sqrt(2*c.massZRelative*9.1093837015e-31*kT*1.602176634e-19)/1.054571817e-34;
    values{s} = 2*reservoir.degeneracy*kzScale/(pi*reservoir.widthM)*value;
end
source = values{1}; drain = values{2};
end
