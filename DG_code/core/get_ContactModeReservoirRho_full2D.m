function [rhoBoundary, info] = get_ContactModeReservoirRho_full2D(mat, p, sideName, Ef, Vxy)
%GET_CONTACTMODERESERVOIRRHO_FULL2D Source/Drain data from contact modes.
%
% This helper implements the more physical Full-2D reservoir model for the
% rho-basis transport equation. At the selected contact edge it first solves
% the transverse contact eigenproblem
%
%   H_perp phi_n(y) = E_n phi_n(y).
%
% The incoming reservoir density matrix is then assembled as
%
%   rho_c(Y,rho_y,rho_x) =
%       sum_n phi_n(Y+rho_y/2) conj(phi_n(Y-rho_y/2)) F_n(rho_x).
%
% Here F_n(rho_x) is the cosine transform of a parabolic longitudinal Fermi
% distribution for the same transverse mode. This deliberately keeps the
% current Full-2D implementation parabolic; a future non-parabolic model can
% replace only the longitudinal occupation and characteristic classification.
% By default the Fermi level is obtained from charge neutrality using these
% same contact modes and their mode-averaged x/z masses. The externally
% supplied Ef is retained as an optional legacy/fallback value.
%
% Relative-vector order follows the rest of the Full-2D code:
%   rho_x is fastest, then rho_y. The returned array therefore has size
%   [N_rho_x*N_rho_y, N_Y_face_DOF].

if nargin < 5 || isempty(Vxy)
    Vxy = readField(mat, 'V', []);
end

params = getDGParams(mat);
constants = physicalConstants;
temp = readField(mat, 'Temp', 300);
degFactor = readDegeneracy(mat);

[contact, contactInfo] = contactProblem(mat, p, sideName, Vxy);
nModes = contactModeCount(mat, params, contact.nGrid);
[modeEnergy, modeVector, modeNormalization] = solveContactModes( ...
    contact.H, nModes, contact.gridY);

modeMassX = modeAveragedMass(modeVector, contact.mXRelative, ...
    readField(mat, 'me_x_ch', 0.041), constants, contact.gridY);
modeMassZ = modeAveragedMass(modeVector, contact.mZRelative, ...
    readField(mat, 'me_z_ch', readField(mat, 'me_x_ch', 0.041)), ...
    constants, contact.gridY);

nKx = readParam(params, 'full2D_contactMode_Nkx', ...
    readParam(params, 'full2D_reservoir_Nkx', p.relative.NrhoX));
[kx, dkx] = makeKGrid(nKx, p.relative.rhoX.length);
cosX = cos(p.relative.rhoX.cells(:)*kx);

integrateKz = shouldIntegrateKz(mat, params);
[EfUsed, fermiLevelInfo] = resolveContactFermiLevel(mat, p, params, ...
    sideName, Ef, modeEnergy, modeMassX, modeMassZ, temp, degFactor, ...
    constants, integrateKz);
if isempty(EfUsed) || ~isscalar(EfUsed) || ~isfinite(EfUsed)
    [rhoBoundary, info] = zeroReservoir(p, sideName, ...
        'No finite contact-mode reservoir Fermi level could be determined.');
    info.fermiLevel = fermiLevelInfo;
    return
end

kzInfo = makeKzInfo(params, modeMassZ, modeEnergy, EfUsed, temp, ...
    constants, integrateKz);
[quadratureDensity, quadratureDensityByMode] = contactQuadratureDensity( ...
    kx, dkx, modeEnergy, modeMassX, modeMassZ, EfUsed, temp, degFactor, ...
    constants, kzInfo);
fermiLevelInfo.quadratureSheetDensity = quadratureDensity;
fermiLevelInfo.quadratureSheetDensityByMode = quadratureDensityByMode;
if integrateKz && isfinite(fermiLevelInfo.targetSheetDensityM2) ...
        && fermiLevelInfo.targetSheetDensityM2 > 0
    fermiLevelInfo.quadratureRelativeResidual = ...
        abs(quadratureDensity-fermiLevelInfo.targetSheetDensityM2) ...
        /fermiLevelInfo.targetSheetDensityM2;
else
    fermiLevelInfo.quadratureRelativeResidual = NaN;
end

yFace = p.dg.Y.nodes(:);
rhoBoundary = zeros(p.relative.nDof, numel(yFace));
rhoXModes = zeros(p.relative.NrhoX, nModes);
transverseNormByMode = zeros(nModes, 1);

for im = 1:nModes
    rhoXMode = transformLongitudinalMode(cosX, kx, dkx, modeEnergy(im), ...
        modeMassX(im), modeMassZ(im), EfUsed, temp, degFactor, constants, kzInfo);
    rhoYMode = transverseModeDensity(modeVector(:, im), ...
        contact.gridY, yFace, p.relative.rhoY.cells(:));

    rhoXModes(:, im) = rhoXMode;
    transverseNormByMode(im) = norm(rhoYMode, 'fro');

    for iy = 1:numel(yFace)
        rhoGrid = rhoXMode*rhoYMode(:, iy).';
        rhoBoundary(:, iy) = rhoBoundary(:, iy) + rhoGrid(:);
    end
end

info = struct;
info.origin = 'contact-modes';
info.model = 'contact-modes';
info.side = sideName;
info.basis = 'rho_x/rho_y';
info.relativeOrder = {'rho_x', 'rho_y'};
info.modeCount = nModes;
info.modeEnergy = modeEnergy(:);
info.modeMassXKg = modeMassX(:);
info.modeMassZKg = modeMassZ(:);
info.modeNormalization = modeNormalization;
info.transverseModeNorm = transverseNormByMode;
info.contact = contactInfo;
info.kxRange = [min(kx), max(kx)];
info.dkx = dkx;
info.kz = kzInfo.info;
info.temperature = temp;
info.degeneracy = degFactor;
info.fermiLevelInput = scalarOrNaN(Ef);
info.fermiLevelUsed = EfUsed;
info.fermiLevel = fermiLevelInfo;
info.size = size(rhoBoundary);
info.note = ['Source/Drain data are built from transverse contact modes ', ...
    'phi_n(y) and their density matrix phi_n(Y+rho_y/2)phi_n^*(Y-rho_y/2). ', ...
    'The longitudinal occupation is parabolic in k_x.'];
info.longitudinalTransform = ...
    'F_n(rho_x) = int f_n(k_x) cos(rho_x*k_x) dk_x/(2*pi).';
info.rhoXModeProfiles = rhoXModes;
end

function [contact, info] = contactProblem(mat, p, sideName, Vxy)
sideIndex = sideToXIndex(mat, sideName);
gridY = getCoordinateVector(mat, 'y', 'dy', 'Ny')*readField(p, 'coordinateScale', 1);
[gridY, sortIds] = sort(gridY(:));

Vfield = fieldToXY(Vxy, mat, 0);
mYfield = fieldToXY(readField(mat, 'me_y', []), mat, readField(mat, 'me_y_ch', 0.041));
mXfield = fieldToXY(readField(mat, 'me_x', []), mat, readField(mat, 'me_x_ch', 0.041));
mZfield = fieldToXY(readField(mat, 'me_z', []), mat, readField(mat, 'me_z_ch', readField(mat, 'me_x_ch', 0.041)));

Vedge = Vfield(sideIndex, sortIds).';
mYRelative = sanitizeRelativeMass(mYfield(sideIndex, sortIds).', readField(mat, 'me_y_ch', 0.041));
mXRelative = sanitizeRelativeMass(mXfield(sideIndex, sortIds).', readField(mat, 'me_x_ch', 0.041));
mZRelative = sanitizeRelativeMass(mZfield(sideIndex, sortIds).', readField(mat, 'me_z_ch', readField(mat, 'me_x_ch', 0.041)));

H = transverseHamiltonian(gridY, Vedge, mYRelative);

contact = struct;
contact.H = H;
contact.gridY = gridY;
contact.Vedge = Vedge;
contact.mYRelative = mYRelative;
contact.mXRelative = mXRelative;
contact.mZRelative = mZRelative;
contact.nGrid = numel(gridY);

info = struct;
info.side = sideName;
info.xIndex = sideIndex;
info.yRange = [gridY(1), gridY(end)];
info.VRange = [min(Vedge), max(Vedge)];
info.meXRelativeRange = [min(mXRelative), max(mXRelative)];
info.meYRelativeRange = [min(mYRelative), max(mYRelative)];
info.meZRelativeRange = [min(mZRelative), max(mZRelative)];
info.HSize = size(H);
info.HNnz = nnz(H);
info.note = ['The transverse contact Hamiltonian is built on the material ', ...
    'Y grid at the selected Source/Drain X edge.'];
end

function H = transverseHamiltonian(gridY, Vedge, mYRelative)
nY = numel(gridY);
if nY < 3
    error('DG:Full2D:ContactModeGridTooSmall', ...
        'Contact-mode reservoir needs at least three Y grid points.');
end

constants = physicalConstants;
dy = mean(diff(gridY));
t0 = constants.hbar^2/(2*constants.m0*constants.q*dy^2);

interior = (2:nY-1).';
leftValue = -2*t0./(mYRelative(interior-1) + mYRelative(interior));
rightValue = -2*t0./(mYRelative(interior+1) + mYRelative(interior));
centerValue = -leftValue - rightValue + Vedge(interior);

% Match the boundary convention used by eval_transversal_Hamiltonian.m so
% the contact-mode model is comparable to the existing Mode-Space setup.
ii = [interior; interior; interior; 1; 1; nY; nY];
jj = [interior-1; interior; interior+1; 1; 2; nY; nY-1];
vv = [leftValue; centerValue; rightValue; ...
    centerValue(1); rightValue(1); centerValue(end); leftValue(end)];
H = sparse(ii, jj, vv, nY, nY);
end

function [modeEnergy, modeVector, normInfo] = solveContactModes(H, nModes, gridY)
nGrid = size(H, 1);
if nModes >= nGrid
    [V, E] = eig(full(H));
else
    try
        [V, E] = eigs(H, nModes, 'sm');
    catch
        [V, E] = eig(full(H));
    end
end

modeEnergy = real(diag(E));
[modeEnergy, ids] = sort(modeEnergy, 'ascend');
modeVector = real(V(:, ids(1:nModes)));
modeEnergy = modeEnergy(1:nModes);

weights = quadratureWeights(gridY);
rawIntegral = zeros(nModes, 1);
normalizedIntegral = zeros(nModes, 1);
for im = 1:nModes
    % EIGS returns vectors with Euclidean norm one. For a density matrix in
    % physical Y coordinates we need continuum normalization,
    % int |phi_n(y)|^2 dy = 1. Without this dy-weighted normalization the
    % contact reservoir scales with the transverse grid spacing.
    rawIntegral(im) = sum(weights.*abs(modeVector(:, im)).^2);
    modeNorm = sqrt(rawIntegral(im));
    if modeNorm > 0
        modeVector(:, im) = modeVector(:, im)/modeNorm;
    end
    if modeVector(end, im) < 0
        modeVector(:, im) = -modeVector(:, im);
    end
    normalizedIntegral(im) = sum(weights.*abs(modeVector(:, im)).^2);
end

normInfo = struct;
normInfo.method = 'trapezoidal-y-integral';
normInfo.rawIntegral = rawIntegral;
normInfo.normalizedIntegral = normalizedIntegral;
normInfo.note = ...
    'Contact modes are normalized with int |phi(y)|^2 dy = 1 to avoid dy-dependent reservoir amplitudes.';
end

function rhoXMode = transformLongitudinalMode(cosX, kx, dkx, modeEnergy, ...
        massX, massZ, Ef, temp, degFactor, constants, kzInfo)
fKx = longitudinalOccupation(kx, modeEnergy, massX, massZ, Ef, temp, ...
    degFactor, constants, kzInfo);
rhoXMode = cosX*fKx*dkx/(2*pi);
end

function fKx = longitudinalOccupation(kx, modeEnergy, massX, massZ, Ef, ...
        temp, degFactor, constants, kzInfo)
EkX = constants.hbar^2*kx(:).^2/(2*massX*constants.q);

if kzInfo.integrate
    kz = kzInfo.k(:).';
    EkZ = constants.hbar^2*kz.^2/(2*massZ*constants.q);
    energy = modeEnergy + EkX + EkZ;
    arg = constants.q*(energy - Ef)/(constants.kB*temp);
    arg = min(max(arg, -80), 80);
    fermi = 1./(1 + exp(arg));
    fKx = 2*degFactor*trapz(kz, fermi, 2)/pi;
else
    energy = modeEnergy + EkX;
    arg = constants.q*(energy - Ef)/(constants.kB*temp);
    arg = min(max(arg, -80), 80);
    fKx = 2*degFactor./(1 + exp(arg));
end
end

function [density, densityByMode] = contactQuadratureDensity(kx, dkx, ...
        modeEnergy, modeMassX, modeMassZ, Ef, temp, degFactor, constants, kzInfo)
nModes = numel(modeEnergy);
densityByMode = zeros(nModes, 1);
for im = 1:nModes
    fKx = longitudinalOccupation(kx, modeEnergy(im), modeMassX(im), ...
        modeMassZ(im), Ef, temp, degFactor, constants, kzInfo);
    densityByMode(im) = sum(real(fKx))*dkx/(2*pi);
end
density = sum(densityByMode);
end

function rhoYMode = transverseModeDensity(modeVector, gridY, yFace, rhoY)
rhoYMode = zeros(numel(rhoY), numel(yFace));
for iy = 1:numel(yFace)
    yPlus = yFace(iy) + 0.5*rhoY;
    yMinus = yFace(iy) - 0.5*rhoY;
    phiPlus = interp1(gridY, modeVector, yPlus, 'linear', 0);
    phiMinus = interp1(gridY, modeVector, yMinus, 'linear', 0);
    rhoYMode(:, iy) = phiPlus.*conj(phiMinus);
end
end

function massKg = modeAveragedMass(modeVector, massRelative, defaultRelative, constants, gridY)
massRelative = sanitizeRelativeMass(massRelative, defaultRelative);
weights = quadratureWeights(gridY);
weight = abs(modeVector).^2;
weightSum = sum(weight.*weights, 1).';
relativeMass = (weight.'*(weights.*massRelative(:)))./max(weightSum, eps);
massKg = sanitizeRelativeMass(relativeMass, defaultRelative)*constants.m0;
end

function weights = quadratureWeights(gridY)
gridY = gridY(:);
if numel(gridY) < 2
    weights = 1;
    return
end

dyLeft = [0; diff(gridY)];
dyRight = [diff(gridY); 0];
weights = 0.5*(dyLeft + dyRight);
weights(1) = 0.5*(gridY(2)-gridY(1));
weights(end) = 0.5*(gridY(end)-gridY(end-1));
end

function nModes = contactModeCount(mat, params, nGrid)
defaultModes = min(4, nGrid-1);
if isfield(mat, 'n_of_modes') && ~isempty(mat.n_of_modes)
    defaultModes = mat.n_of_modes;
end
nModes = round(readParam(params, 'full2D_contactModeCount', defaultModes));
nModes = max(1, min(nModes, nGrid));
end

function [EfUsed, info] = resolveContactFermiLevel(mat, p, params, sideName, ...
        EfInput, modeEnergy, modeMassX, modeMassZ, temp, degFactor, ...
        constants, integrateKz)
requestedModel = normalizeFermiModel(readParam(params, ...
    'full2D_contactModeFermiModel', 'contact-neutrality'));
[targetDensity, targetInfo] = contactTargetSheetDensity( ...
    mat, p, params, sideName);
inputValue = scalarOrNaN(EfInput);

info = struct;
info.requestedModel = requestedModel;
info.model = requestedModel;
info.inputEV = inputValue;
info.usedEV = NaN;
info.targetSheetDensityM2 = targetDensity;
info.target = targetInfo;
info.achievedSheetDensityM2 = NaN;
info.achievedSheetDensityByModeM2 = NaN(numel(modeEnergy), 1);
info.relativeResidual = NaN;
info.iterations = 0;
info.converged = false;
info.note = '';

if strcmp(requestedModel, 'contact-neutrality')
    if ~integrateKz
        info.model = 'external-fallback';
        info.note = ['Contact neutrality against a volumetric doping target ', ...
            'requires the unresolved k_z continuum. Using the supplied Ef.'];
    elseif ~isfinite(targetDensity) || targetDensity <= 0
        info.model = 'external-fallback';
        info.note = ['No positive contact sheet-density target is available. ', ...
            'Using the supplied Ef.'];
    else
        [EfUsed, solveInfo] = solveContactNeutrality(modeEnergy, ...
            modeMassX, modeMassZ, targetDensity, temp, degFactor, ...
            constants, params, inputValue);
        info.usedEV = EfUsed;
        info.achievedSheetDensityM2 = solveInfo.achievedDensity;
        info.achievedSheetDensityByModeM2 = solveInfo.densityByMode;
        info.relativeResidual = solveInfo.relativeResidual;
        info.iterations = solveInfo.iterations;
        info.converged = solveInfo.converged;
        info.bracketEV = solveInfo.bracketEV;
        info.note = ['Ef was solved from the same transverse contact modes ', ...
            'and x/z density-of-states masses used by the reservoir.'];
        return
    end
end

EfUsed = inputValue;
info.usedEV = EfUsed;
if isfinite(EfUsed) && integrateKz
    [density, densityByMode] = contactSheetDensity(EfUsed, modeEnergy, ...
        modeMassX, modeMassZ, temp, degFactor, constants);
    info.achievedSheetDensityM2 = density;
    info.achievedSheetDensityByModeM2 = densityByMode;
    if isfinite(targetDensity) && targetDensity > 0
        info.relativeResidual = abs(density-targetDensity)/targetDensity;
    end
end
end

function model = normalizeFermiModel(rawModel)
model = lower(strrep(strrep(strtrim(char(rawModel)), '_', '-'), ' ', '-'));
switch model
    case {'contact', 'contact-modes', 'neutrality', 'contact-neutrality', 'solve'}
        model = 'contact-neutrality';
    case {'external', 'legacy', 'supplied', 'input'}
        model = 'external';
    otherwise
        error('DG:Full2D:UnknownContactFermiModel', ...
            ['Unknown full2D_contactModeFermiModel "%s". Use ', ...
             '"contact-neutrality" or "external".'], char(rawModel));
end
end

function [targetDensity, info] = contactTargetSheetDensity( ...
        mat, p, params, sideName)
sideName = lower(char(sideName));
overrideName = ['full2D_', sideName, 'ContactSheetDensity'];
targetDensity = readParam(params, overrideName, NaN);

info = struct;
info.side = sideName;
info.units = 'm^-2';
info.origin = overrideName;
info.dopingCm3 = NaN;
info.contactWidthM = NaN;

if isfiniteScalar(targetDensity)
    validateattributes(targetDensity, {'numeric'}, ...
        {'real', 'finite', 'scalar', 'positive'}, mfilename, overrideName);
    return
end

switch sideName
    case 'source'
        dopingName = 'N_s';
    case 'drain'
        dopingName = 'N_d';
    otherwise
        targetDensity = NaN;
        info.origin = 'unavailable';
        return
end

widthName = ['full2D_', sideName, 'ContactWidth'];
if isfield(mat, dopingName) && ~isempty(mat.(dopingName)) ...
        && isfield(mat, 'W_c') && ~isempty(mat.W_c)
    doping = mat.(dopingName);
    widthGridUnits = readParam(params, widthName, mat.W_c);
    validateattributes(doping, {'numeric'}, ...
        {'real', 'finite', 'scalar', 'positive'}, mfilename, dopingName);
    validateattributes(widthGridUnits, {'numeric'}, ...
        {'real', 'finite', 'scalar', 'positive'}, mfilename, widthName);
    widthM = widthGridUnits*readField(p, 'coordinateScale', 1);
    targetDensity = doping*1e6*widthM;
    info.origin = [dopingName, '*', widthName];
    info.dopingCm3 = doping;
    info.contactWidthM = widthM;
else
    targetDensity = NaN;
    info.origin = 'unavailable';
end
end

function [Ef, info] = solveContactNeutrality(modeEnergy, modeMassX, ...
        modeMassZ, targetDensity, temp, degFactor, constants, params, EfGuess)
relativeTolerance = readParam(params, ...
    'full2D_contactModeFermiTolerance', 1e-10);
maxIterations = readParam(params, ...
    'full2D_contactModeFermiMaxIterations', 200);
validateattributes(relativeTolerance, {'numeric'}, ...
    {'real', 'finite', 'scalar', 'positive'}, mfilename, ...
    'full2D_contactModeFermiTolerance');
validateattributes(maxIterations, {'numeric'}, ...
    {'real', 'finite', 'scalar', 'integer', 'positive'}, mfilename, ...
    'full2D_contactModeFermiMaxIterations');

kTEV = constants.kB*temp/constants.q;
modeMinimum = min(real(modeEnergy(:)));
step = max(0.5, 20*kTEV);
lower = modeMinimum - max(2, 80*kTEV);
upper = modeMinimum + max(2, 80*kTEV);
if isfinite(EfGuess)
    lower = min(lower, EfGuess-step);
    upper = max(upper, EfGuess+step);
end

[densityLower, ~] = contactSheetDensity(lower, modeEnergy, modeMassX, ...
    modeMassZ, temp, degFactor, constants);
[densityUpper, ~] = contactSheetDensity(upper, modeEnergy, modeMassX, ...
    modeMassZ, temp, degFactor, constants);
for iexpand = 1:50
    if densityLower <= targetDensity && densityUpper >= targetDensity
        break
    end
    step = 2*step;
    if densityLower > targetDensity
        lower = lower-step;
        [densityLower, ~] = contactSheetDensity(lower, modeEnergy, ...
            modeMassX, modeMassZ, temp, degFactor, constants);
    end
    if densityUpper < targetDensity
        upper = upper+step;
        [densityUpper, ~] = contactSheetDensity(upper, modeEnergy, ...
            modeMassX, modeMassZ, temp, degFactor, constants);
    end
end
if densityLower > targetDensity || densityUpper < targetDensity
    error('DG:Full2D:ContactFermiBracketFailed', ...
        'Could not bracket the contact-mode Fermi level for %g m^-2.', ...
        targetDensity);
end

converged = false;
iterations = 0;
for iterations = 1:maxIterations
    midpoint = 0.5*(lower+upper);
    [densityMidpoint, ~] = contactSheetDensity(midpoint, modeEnergy, ...
        modeMassX, modeMassZ, temp, degFactor, constants);
    relativeResidual = abs(densityMidpoint-targetDensity)/targetDensity;
    if relativeResidual <= relativeTolerance
        lower = midpoint;
        upper = midpoint;
        converged = true;
        break
    elseif densityMidpoint < targetDensity
        lower = midpoint;
    else
        upper = midpoint;
    end
end

Ef = 0.5*(lower+upper);
[achievedDensity, densityByMode] = contactSheetDensity(Ef, modeEnergy, ...
    modeMassX, modeMassZ, temp, degFactor, constants);
info = struct;
info.achievedDensity = achievedDensity;
info.densityByMode = densityByMode;
info.relativeResidual = abs(achievedDensity-targetDensity)/targetDensity;
info.iterations = iterations;
info.converged = converged || info.relativeResidual <= relativeTolerance;
info.bracketEV = [lower, upper];
end

function [density, densityByMode] = contactSheetDensity(Ef, modeEnergy, ...
        modeMassX, modeMassZ, temp, degFactor, constants)
kBT = constants.kB*temp;
eta = constants.q*(Ef-real(modeEnergy(:)))/kBT;
dosMass = sqrt(real(modeMassX(:)).*real(modeMassZ(:)));
prefactor = degFactor*dosMass*kBT/(pi*constants.hbar^2);
densityByMode = prefactor.*stableSoftplus(eta);
density = sum(densityByMode);
end

function value = stableSoftplus(argument)
value = zeros(size(argument));
large = argument > 40;
small = argument < -40;
middle = ~(large | small);
value(large) = argument(large);
value(small) = exp(argument(small));
value(middle) = log1p(exp(argument(middle)));
end

function value = scalarOrNaN(rawValue)
value = NaN;
if isnumeric(rawValue) && isscalar(rawValue) && isfinite(rawValue)
    value = real(rawValue);
end
end

function tf = isfiniteScalar(value)
tf = isnumeric(value) && isscalar(value) && isfinite(value);
end

function tf = shouldIntegrateKz(mat, params)
defaultValue = isfield(mat, 'me_z') || isfield(mat, 'me_z_ch');
tf = readLogicalParam(params, 'full2D_contactModeIntegrateKz', defaultValue);
end

function kzInfo = makeKzInfo(params, modeMassZ, modeEnergy, Ef, temp, ...
        constants, integrateKz)
kzInfo = struct;
kzInfo.integrate = integrateKz;
kzInfo.k = 0;
kzInfo.info = struct;
kzInfo.info.integrate = integrateKz;

if ~integrateKz
    kzInfo.info.note = 'No k_z integration; occupation is a 1D parabolic Fermi factor in k_x.';
    return
end

nKz = readParam(params, 'full2D_contactMode_Nkz', ...
    readParam(params, 'full2D_reservoir_Nkz', 1001));
minimumEnergyWindow = readParam(params, ...
    'full2D_contactModeKzEnergyWindow', 0.3);
thermalTailKBT = readParam(params, ...
    'full2D_contactModeKzThermalTail', 10);

validateattributes(minimumEnergyWindow, {'numeric'}, ...
    {'real', 'finite', 'scalar', 'nonnegative'}, mfilename, ...
    'full2D_contactModeKzEnergyWindow');
validateattributes(thermalTailKBT, {'numeric'}, ...
    {'real', 'finite', 'scalar', 'nonnegative'}, mfilename, ...
    'full2D_contactModeKzThermalTail');

% The transverse modes replace the explicit k_y integration, but k_z is
% still a continuous reservoir coordinate. Its kinetic-energy interval must
% extend beyond the occupied part of the lowest contact mode. A fixed cutoff
% can otherwise truncate the Fermi sea when Ef-E_mode is comparable to the
% configured window. Keep the user value as a lower bound and add a thermal
% tail so the occupation at the numerical endpoint is exponentially small.
thermalTailEV = thermalTailKBT*constants.kB*temp/constants.q;
occupiedWindow = max(max(real(Ef(:))) - min(real(modeEnergy(:))), 0);
adaptiveEnergyWindow = occupiedWindow + thermalTailEV;
energyWindow = max(minimumEnergyWindow, adaptiveEnergyWindow);

nKz = max(2, round(nKz));
mzRef = max(real(modeMassZ(:)));
kzMax = sqrt(max(energyWindow, eps)*2*mzRef*constants.q)/constants.hbar;
kzInfo.k = linspace(0, kzMax, nKz);
kzInfo.info.nKz = nKz;
kzInfo.info.energyWindowEV = energyWindow;
kzInfo.info.minimumEnergyWindowEV = minimumEnergyWindow;
kzInfo.info.adaptiveEnergyWindowEV = adaptiveEnergyWindow;
kzInfo.info.occupiedEnergyWindowEV = occupiedWindow;
kzInfo.info.thermalTailEV = thermalTailEV;
kzInfo.info.thermalTailKBT = thermalTailKBT;
kzInfo.info.range = [0, kzMax];
kzInfo.info.note = ['k_z is integrated assuming parabolic dispersion and ', ...
    'translational invariance in the unmodeled width direction. The ', ...
    'kinetic-energy cutoff is max(user minimum, Ef-min(E_mode)+thermal tail).'];
end

function [k, dk] = makeKGrid(nK, lengthValue)
nK = max(2, round(nK));
k = -2*pi/lengthValue*((1:nK) - 0.5*(nK+1));
dk = abs(k(2)-k(1));
end

function field = fieldToXY(rawField, mat, defaultValue)
if isempty(rawField)
    field = defaultValue*ones(mat.Nx, mat.Ny);
    return
end

rawField = squeeze(rawField);
fieldSize = size(rawField);

if ismatrix(rawField) && isequal(fieldSize, [mat.Nx, mat.Ny])
    field = rawField;
elseif ismatrix(rawField) && isequal(fieldSize, [mat.Ny, mat.Nx])
    field = rawField.';
elseif isvector(rawField) && numel(rawField) == mat.Nx
    field = repmat(rawField(:), 1, mat.Ny);
elseif isvector(rawField) && numel(rawField) == mat.Ny
    field = repmat(rawField(:).', mat.Nx, 1);
elseif ndims(rawField) == 3 && size(rawField, 2) == mat.Nx && size(rawField, 3) == mat.Ny
    field = squeeze(rawField(1, :, :));
elseif ndims(rawField) == 4 && size(rawField, 2) == mat.Nx && size(rawField, 3) == mat.Ny
    zIndex = ceil(size(rawField, 4)/2);
    field = squeeze(rawField(1, :, :, zIndex));
else
    error('DG:Full2D:InvalidMaterialField', ...
        'Cannot interpret material field with size [%s] as an X/Y field.', ...
        num2str(size(rawField)));
end
end

function sideIndex = sideToXIndex(mat, sideName)
switch lower(sideName)
    case 'source'
        sideIndex = 1;
    case 'drain'
        sideIndex = mat.Nx;
    otherwise
        error('DG:Full2D:UnknownReservoirSide', ...
            'Unknown reservoir side "%s". Use "source" or "drain".', sideName);
end
end

function coord = getCoordinateVector(mat, coordName, spacingName, countName)
if isfield(mat, coordName) && ~isempty(mat.(coordName))
    coord = mat.(coordName)(:).';
else
    coord = (0:mat.(countName)-1)*mat.(spacingName);
end
end

function mass = sanitizeRelativeMass(mass, defaultValue)
mass = real(mass(:));
mass(~isfinite(mass) | mass <= 0) = defaultValue;
end

function params = getDGParams(mat)
params = struct;
if isfield(mat, 'dg') && isfield(mat.dg, 'params')
    params = mat.dg.params;
end
end

function value = readParam(params, name, defaultValue)
value = defaultValue;
if isstruct(params)
    if isfield(params, name) && ~isempty(params.(name))
        value = params.(name);
    elseif isfield(params, 'full2D') && isstruct(params.full2D) ...
            && isfield(params.full2D, name) && ~isempty(params.full2D.(name))
        value = params.full2D.(name);
    end
end
end

function value = readField(s, name, defaultValue)
value = defaultValue;
if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
    value = s.(name);
end
end

function value = readLogicalParam(params, name, defaultValue)
rawValue = readParam(params, name, defaultValue);
if islogical(rawValue)
    value = rawValue;
elseif isnumeric(rawValue)
    value = rawValue ~= 0;
else
    value = any(strcmpi(char(rawValue), {'true', 'yes', 'on', '1'}));
end
end

function degFactor = readDegeneracy(mat)
degFactor = 1;
if isfield(mat, 'deg_factor') && ~isempty(mat.deg_factor)
    degFactor = mat.deg_factor(1);
end
end

function c = physicalConstants
c.hbar = 1.054571817e-34;
c.q = 1.602176634e-19;
c.kB = 1.38064852e-23;
c.m0 = 9.1093837015e-31;
end

function [rhoBoundary, info] = zeroReservoir(p, sideName, note)
rhoBoundary = zeros(p.relative.nDof, p.dg.Y.nDof);
info = struct;
info.origin = 'zero-placeholder';
info.model = 'contact-modes';
info.side = sideName;
info.note = note;
info.size = size(rhoBoundary);
end
