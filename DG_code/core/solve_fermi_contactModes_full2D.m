function [EfL, EfR, info] = solve_fermi_contactModes_full2D( ...
        mat, Vxy, EfGuessL, EfGuessR)
%SOLVE_FERMICONTACTMODES_FULL2D Equilibrium contact Ef from Full-2D modes.
%
% The old Mode-Space Fermi levels may be supplied as numerical guesses. The
% returned values are the levels actually selected by the contact-mode
% reservoir model, i.e. from charge neutrality with the newly calculated
% transverse Source/Drain modes. This helper is used before the Poisson loop
% so the gate work-function boundary and the transport reservoirs share the
% same equilibrium energy reference.

if nargin < 4
    EfGuessR = [];
end
if nargin < 3
    EfGuessL = [];
end
if nargin < 2 || isempty(Vxy)
    Vxy = mat.V;
end

p = initParams_full2D(mat, Vxy, EfGuessL, EfGuessR);
[~, sourceInfo] = get_ContactModeReservoirRho_full2D( ...
    mat, p, 'source', EfGuessL, Vxy);
[~, drainInfo] = get_ContactModeReservoirRho_full2D( ...
    mat, p, 'drain', EfGuessR, Vxy);

EfL = selectedFermiLevel(sourceInfo, EfGuessL, 'source');
EfR = selectedFermiLevel(drainInfo, EfGuessR, 'drain');

info = struct;
info.source = sourceInfo.fermiLevel;
info.drain = drainInfo.fermiLevel;
info.EfL = EfL;
info.EfR = EfR;
info.note = ['Equilibrium Full-2D contact levels evaluated before the ', ...
    'self-consistent loop from the contact-mode neutrality model.'];
end

function Ef = selectedFermiLevel(reservoirInfo, fallback, sideName)
Ef = fallback;
if isfield(reservoirInfo, 'fermiLevelUsed') ...
        && isscalar(reservoirInfo.fermiLevelUsed) ...
        && isfinite(reservoirInfo.fermiLevelUsed)
    Ef = reservoirInfo.fermiLevelUsed;
end
if isempty(Ef) || ~isscalar(Ef) || ~isfinite(Ef)
    error('DG:Full2D:ContactFermiUnavailable', ...
        'No finite Full-2D contact Fermi level is available for %s.', sideName);
end
end
