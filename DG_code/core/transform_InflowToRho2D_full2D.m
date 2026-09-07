function [rho2D, info] = transform_InflowToRho2D_full2D(inflow1D, p, opts)
%TRANSFORM_INFLOWTORHO2D_FULL2D Lift existing 1D inflow data to rho_x/rho_y.
%
% If opts.inputBasis is 'rho_x', inflow1D is already a longitudinal rho_x
% vector. If it is 'phase_x', p.relative.PhiX must be supplied and the data
% are transformed by rho_x = PhiX*c_x, matching the old 1D phase-to-rho idea.
%
% The transverse rho_y profile is separated deliberately. If no profile is
% provided, a discrete even delta profile at rho_y=0 is used as a neutral
% placeholder until a true 2D reservoir distribution in (k_x,k_y) is added.

if nargin < 3 || isempty(opts)
    opts = struct;
end
if ~isfield(opts, 'inputBasis')
    opts.inputBasis = 'rho_x';
end

inflow1D = inflow1D(:);
switch lower(opts.inputBasis)
    case 'rho_x'
        rhoX = inflow1D;
    case 'phase_x'
        if ~isfield(p.relative, 'PhiX')
            error('DG:Full2D:MissingPhiX', ...
                'p.relative.PhiX is required to transform phase_x inflow data.');
        end
        rhoX = p.relative.PhiX*inflow1D;
    otherwise
        error('DG:Full2D:UnknownInflowBasis', ...
            'Unknown inflow basis "%s". Use "rho_x" or "phase_x".', opts.inputBasis);
end

if numel(rhoX) ~= p.relative.NrhoX
    error('DG:Full2D:InvalidInflowLength', ...
        'Longitudinal inflow has %d entries, expected %d rho_x cells.', ...
        numel(rhoX), p.relative.NrhoX);
end

if isfield(opts, 'rhoYProfile') && ~isempty(opts.rhoYProfile)
    rhoYProfile = opts.rhoYProfile(:);
else
    rhoYProfile = zeros(p.relative.NrhoY, 1);
    [~, centerId] = min(abs(p.relative.rhoY.cells));
    rhoYProfile(centerId) = 1;
end

if numel(rhoYProfile) ~= p.relative.NrhoY
    error('DG:Full2D:InvalidTransverseProfile', ...
        'rhoYProfile has %d entries, expected %d rho_y cells.', ...
        numel(rhoYProfile), p.relative.NrhoY);
end

% Relative-vector order is (rho_x, rho_y), with rho_x fastest.
rho2D = kron(rhoYProfile, rhoX);

info = struct;
info.inputBasis = opts.inputBasis;
info.outputBasis = 'rho_x/rho_y';
info.rhoYProfile = rhoYProfile;
info.relativeOrder = {'rho_x', 'rho_y'};
end
