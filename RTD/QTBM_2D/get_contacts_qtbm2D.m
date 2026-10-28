function contacts = get_contacts_qtbm2D(mat,potentialEV,options)
%GET_CONTACTS_QTBM2D Read-only adapter to the existing periodic mode helper.
% Restricted to Y-homogeneous fields and globally constant effective masses.
% energyEV includes the local contact band edge. The transverse offsets are
% stored separately to avoid counting that edge twice in the longitudinal H.
validateattributes(potentialEV,{'double','single'}, ...
    {'real','finite','size',[mat.Nx,mat.Ny]});
if max(abs(potentialEV-potentialEV(:,1)),[],'all') > 1e-10
    error('QTBM2D:NonseparablePotential','Separable transport requires V(X,Y)=V(X).');
end
fields = {'Vc','eps','Nd','N_d','N_a'};
for k = 1:numel(fields)
    if isfield(mat,fields{k})
        value = mat.(fields{k});
        if max(abs(value-value(:,1)),[],'all') > 1e-12*max(1,max(abs(value(:))))
            error('QTBM2D:NonseparableDevice','mat.%s must be independent of Y.',fields{k});
        end
    end
end
rtdDir = fileparts(fileparts(mfilename('fullpath')));
oldPath = path;
restorePath = onCleanup(@() path(oldPath));
addpath(fullfile(fileparts(rtdDir),'DG_code','core'));

% Only this local compatibility copy knows the legacy helper's field names.
% NrhoY is a validation allowance, NOT a QTBM relative-coordinate grid.
modeMat = mat;
modeMat.dg.params = struct('full2D_periodicModeOrder',options.modeOrder);
p.coordinateScale = 1e-9;
p.domain.Y = [mat.y(1),mat.y(end)]*p.coordinateScale;
p.relative.NrhoY = max(4,4*options.modeOrder+2);
contacts = struct('widthM',diff(p.domain.Y),'modeOrder',options.modeOrder);
for sideCell = {'source','drain'}
    side = sideCell{1};
    [contact,info,energy,modes,normalization] = ...
        get_PeriodicContactModes_full2D(modeMat,p,side,potentialEV);
    contacts.(side) = struct('energyEV',energy, ...
        'transverseEnergyEV',energy-contact.Vedge(1), ...
        'bandEdgeEV',contact.Vedge(1),'modes',modes, ...
        'massXRelative',contact.mXRelative(1), ...
        'massYRelative',contact.mYRelative(1), ...
        'massZRelative',contact.mZRelative(1), ...
        'normalization',normalization,'info',info);
end
contacts.occupationStatus = 'neutrality-computed-by-get_occupation_qtbm2D';
end
