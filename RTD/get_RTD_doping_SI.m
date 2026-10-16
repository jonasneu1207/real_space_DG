function [Nd,Na] = get_RTD_doping_SI(mat)
%GET_RTD_DOPING_SI Fully ionized donors/acceptors in m^-3, with unit checks.
% Nd/Na are SI; legacy N_d/N_a are cm^-3. If both exist, they must agree.
arguments
    mat (1,1) struct
end
Nd = densityField(mat,'Nd','N_d',true);
Na = densityField(mat,'Na','N_a',false);
end

function value = densityField(mat,siName,legacyName,required)
hasSI = isfield(mat,siName);
hasLegacy = isfield(mat,legacyName);
if ~hasSI && ~hasLegacy
    if required, error('RTD:MissingDoping','Provide mat.Nd [m^-3] or mat.N_d [cm^-3].'); end
    value = zeros(mat.Nx,mat.Ny);
    return
end
if hasSI
    value = mat.(siName);
    validateattributes(value,{'double','single'},{'real','finite','nonnegative','size',[mat.Nx,mat.Ny]});
else
    value = mat.(legacyName)*1e6;
end
if hasLegacy
    legacy = mat.(legacyName);
    validateattributes(legacy,{'double','single'},{'real','finite','nonnegative','size',[mat.Nx,mat.Ny]});
    if hasSI && max(abs(value(:)-legacy(:)*1e6)) > 1e-12*max([1;value(:)])
        error('RTD:InconsistentDopingUnits','mat.%s and mat.%s disagree after cm^-3 to m^-3 conversion.',siName,legacyName);
    end
end
value = double(value);
end
