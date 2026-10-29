function [mat,reference] = prepare_RTD_neumann_contacts(mat)
%PREPARE_RTD_NEUMANN_CONTACTS Calibrate DG modes once, then use fixed mu.
% Flat material bands define the energy reference, NOT the initial phi.
% At bias U: muS=neutralMuSource, muD=muS-U (eV / V convention).
% Unequal neutral contacts imply phiD-phiS=U+builtInVoltage initially.
if ~strcmp(get_YBoundaryType_full2D(mat),'periodic') ...
        || ~strcmpi(mat.dg.params.full2D_reservoirModel,'contact-modes')
    error('RTD:NeumannContactModel','Neumann RTD requires periodic Y and contact-modes reservoirs.');
end
[Nd,Na] = get_RTD_doping_SI(mat);
mat.dg.params.full2D_contactModeFermiModel = 'contact-neutrality';
mat.dg.params.full2D_contactModeIntegrateKz = true;
mat.dg.params.full2D_sourceContactSheetDensity = trapz(mat.y*1e-9,Nd(1,:)-Na(1,:));
mat.dg.params.full2D_drainContactSheetDensity = trapz(mat.y*1e-9,Nd(end,:)-Na(end,:));
p = initParams_full2D(mat,mat.Vc,[],[]);
[~,source] = get_ContactModeReservoirRho_full2D(mat,p,'source',[],mat.Vc);
[~,drain] = get_ContactModeReservoirRho_full2D(mat,p,'drain',[],mat.Vc);
if ~source.fermiLevel.converged || ~drain.fermiLevel.converged
    error('RTD:ReservoirCalibration','Both contact-mode neutrality calibrations must converge.');
end
reference = struct('sourceMuEV',source.fermiLevelUsed, ...
    'neutralDrainMuEV',drain.fermiLevelUsed, ...
    'builtInVoltage',drain.fermiLevelUsed-source.fermiLevelUsed, ...
    'sourceCalibration',source.fermiLevel,'drainCalibration',drain.fermiLevel);
mat.dg.params.full2D_contactModeFermiModel = 'external';
mat.reservoirReference = reference;
end
