function [contact, info, energy, modes, normalization] = ...
        get_PeriodicContactModes_full2D(mat,p,side,Vxy)
%GET_PERIODICCONTACTMODES_FULL2D Fourier-Galerkin modes on a periodic Y cell.
% Includes the complete +/- order pairs. For a homogeneous contact these
% are exact plane waves, including n=0. For a Y-dependent periodic potential
% diagonalize H in the same Fourier basis. Constant mass is required.
% full2D_periodicModeOrder=K retains 2K+1 Fourier basis functions. All their
% eigenstates enter the reservoir; legacy contactModeCount is not used.
params = mat.dg.params;
gridY = mat.y(:)*p.coordinateScale;
width = diff(p.domain.Y);
nUnique = numel(gridY)-1;
dy = width/nUnique;
if nUnique < 2 || max(abs(diff(gridY)-dy)) > 1e-10*dy
    error('DG:Full2D:PeriodicContactGrid', ...
        'Periodic contact modes require a uniform material Y grid including both cell endpoints.');
end
maxOrder = min(floor((nUnique-1)/2),floor((p.relative.NrhoY-2)/4));
order = maxOrder;
if isfield(params,'full2D_periodicModeOrder')
    order = params.full2D_periodicModeOrder;
end
if ~isscalar(order) || ~isfinite(order) || order < 0 || order ~= fix(order) || order > maxOrder
    error('DG:Full2D:PeriodicModeOrder', ...
        'periodicModeOrder must be an integer in 0:%d; refine material Y and/or rho_y for more modes.',maxOrder);
end
ix = 1;
if strcmpi(side,'drain')
    ix = mat.Nx;
end
Vfield = reshape(Vxy,mat.Nx,mat.Ny);
edge = Vfield(ix,:).';
if abs(edge(1)-edge(end)) > 1e-10*max(1,max(abs(edge)))
    error('DG:Full2D:NonperiodicPotential','Contact potential endpoints must agree.');
end
mx = constantMass(mat,'x');
my = constantMass(mat,'y');
mz = constantMass(mat,'z');
hbar = 1.054571817e-34;
q = 1.602176634e-19;
m0 = 9.1093837015e-31;
k = (2*pi/width)*(-order:order);
basis = exp(1i*(gridY(1:end-1)-gridY(1))*k)/sqrt(width);
H = diag(hbar^2*k.^2/(2*my*m0*q)) + dy*(basis'*(edge(1:end-1).*basis));
H = (H+H')/2;
if max(abs(edge-edge(1))) < 1e-13*max(1,abs(edge(1)))
    % Preserve exact +/- plane waves instead of arbitrary degenerate pairs.
    energy = edge(1)+hbar^2*k(:).^2/(2*my*m0*q);
    [energy,ids] = sort(energy);
    coefficients = eye(numel(k));
    coefficients = coefficients(:,ids);
else
    [coefficients,energy] = eig(H,'vector');
    [energy,ids] = sort(real(energy));
    coefficients = coefficients(:,ids);
end
modes = exp(1i*(gridY-gridY(1))*k)*coefficients/sqrt(width);
contact = struct('H',H,'gridY',gridY,'Vedge',edge, ...
    'mXRelative',mx*ones(size(gridY)), ...
    'mYRelative',my*ones(size(gridY)), ...
    'mZRelative',mz*ones(size(gridY)), ...
    'waveNumbers',k,'coefficients',coefficients);
normalization = struct('method','periodic Fourier orthonormal basis', ...
    'rawIntegral',ones(numel(k),1), ...
    'normalizedIntegral',trapz(gridY,abs(modes).^2,1).');
info = struct('side',side,'xIndex',ix,'yRange',p.domain.Y, ...
    'modeExtension','periodic-Fourier-evaluation','period',width, ...
    'fourierOrder',order,'waveNumbers',k,'HSize',size(H), ...
    'note','All retained Fourier states, including both members of every +/- pair, are occupied with their own energies.');
end

function m = constantMass(mat,axis)
m0 = 9.1093837015e-31;
reference = ['me_',axis,'_ch'];
field = ['me_',axis];
if isfield(mat,reference)
    m = mat.(reference);
else
    m = mat.me_x_ch;
end
if m < 1e-25
    m = m/m0;
end
if isfield(mat,field) && ~isempty(mat.(field))
    values = mat.(field);
    values(abs(values)<1e-25) = values(abs(values)<1e-25)/m0;
    if any(~isfinite(values(:))) || max(abs(values(:)-m)) > 1e-10*m
        error('DG:Full2D:PeriodicContactMass', ...
            'Periodic contact modes currently require constant %s matching the reference mass.',field);
    end
end
if ~isscalar(m) || ~isfinite(m) || m <= 0
    error('DG:Full2D:PeriodicContactMass','Periodic contact masses must be positive scalars.');
end
end
