function reservoir = get_occupation_qtbm2D(mat,contacts,energyEV,options)
%GET_OCCUPATION_QTBM2D Neutral periodic-Y contacts, continuous kz and spin 2.
% Neutrality uses the same continuum mx/mz sheet-density formula as DG.
% energyEV is LONGITUDINAL energy including the band edge [eV].
% Optional contacts.fixedMuEV=[muSource,muDrain] bypasses recalibration;
% the current mode energies still enter all occupations and diagnostics.
q = 1.602176634e-19; kB = 1.38064852e-23;
hbar = 1.054571817e-34; m0 = 9.1093837015e-31;
deg = mat.deg_factor; kT = kB*mat.Temp/q;
reservoir = struct('temperatureK',mat.Temp,'degeneracy',deg, ...
    'widthM',contacts.widthM,'kzPoints',options.kzPoints,'kTEV',kT);
fixedMu = isfield(contacts,'fixedMuEV');
if fixedMu
    validateattributes(contacts.fixedMuEV,{'numeric'},{'real','finite','vector','numel',2});
end
for sideCell = {'source','drain'}
    side = sideCell{1}; c = contacts.(side);
    row = 1; if strcmp(side,'drain'), row = mat.Nx; end
    if isfield(mat,'Na'), acceptors = mat.Na(row,:); else, acceptors = mat.N_a(row,:)*1e6; end
    target = trapz(mat.y*1e-9,mat.Nd(row,:)-acceptors);
    if target<=0
        error('QTBM2D:ContactDoping','Electron reservoirs require positive net contact doping.');
    end
    prefactor = deg*sqrt(c.massXRelative*c.massZRelative)*m0*kB*mat.Temp/(pi*hbar^2);
    density = @(mu) prefactor*sum(softplus((mu-c.energyEV(:))/kT));
    if fixedMu
        c.muEV = contacts.fixedMuEV(1+strcmp(side,'drain'));
        c.fermiModel = 'fixed-reservoir';
    else
        lo = min(c.energyEV)-max(2,80*kT); hi = max(c.energyEV)+max(2,80*kT);
        while density(lo)>target, lo = lo-2; end
        while density(hi)<target, hi = hi+2; end
        for iteration = 1:200
            mid = (lo+hi)/2;
            if density(mid)>target, hi = mid; else, lo = mid; end
            if hi-lo<1e-12, break; end
        end
        c.muEV = (lo+hi)/2;
        c.fermiModel = 'contact-neutrality';
    end
    c.targetSheetDensity = target;
    c.achievedSheetDensity = density(c.muEV);
    c.neutralityRelativeError = abs(density(c.muEV)-target)/target;
    reservoir.(side) = c;
end
if ~isempty(energyEV)
    [reservoir.sourceWeightPerArea,reservoir.drainWeightPerArea] = ...
        contact_supply_qtbm2D(energyEV,reservoir);
end
end

function value = softplus(x)
value = max(x,0)+log1p(exp(-abs(x)));
end
