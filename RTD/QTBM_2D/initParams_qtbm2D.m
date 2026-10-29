function options = initParams_qtbm2D(input)
%INITPARAMS_QTBM2D Independent QTBM configuration; see README for units.
if nargin<1, input = struct; end
validateattributes(input,{'struct'},{'scalar'});
options = struct('simulation','self-consistent','Ly',10,'dy',0.25,'bias',(0.0:0.02:0.6), ...
    'modeOrder',2,'energyStepEV',1e-3,'energyRelativeTolerance',2e-4, ...
    'energyTailKT',24,'maxEnergyRefinements',12,'maxEnergyEvaluations',500000, ...
    'kzPoints',64,'currentTolerance',1e-6,'boundStateOccupation','error', ...
    'poissonMethod','newton-raphson','poissonXBoundary','neumann', ...
    'chargeNeutralityTolerance',1e-3,'maxIterations',150,'mixing',1.0, ...
    'andersonDepth',5,'maxPotentialStepV',0.03,'initialPhi',[],'couplingMethod','anderson', ...
    'potentialTolerance',5e-4,'densityTolerance',5e-4, ...
    'poissonResidualTolerance',1,'verbose',true);
names = fieldnames(input);
for k = 1:numel(names)
    name = names{k};
    if ~isfield(options,name), error('QTBM2D:UnknownOption','Unknown QTBM option: %s.',name); end
    options.(name) = input.(name);
end
options.simulation = validatestring(options.simulation,{'prepare','fixed-potential','self-consistent'});
options.poissonMethod = validatestring(options.poissonMethod,{'direct','newton-raphson'});
options.poissonXBoundary = validatestring(options.poissonXBoundary,{'dirichlet','neumann'});
if strcmp(options.poissonXBoundary,'neumann') && strcmp(options.poissonMethod,'direct')
    error('RTD:NeumannNeedsNewton','Neumann X requires newton-raphson.');
end
options.couplingMethod = validatestring(options.couplingMethod,{'anderson','broyden'});
options.boundStateOccupation = validatestring(options.boundStateOccupation,{'error','source','drain','average'});
positive = {'Ly','dy','energyStepEV','energyRelativeTolerance','energyTailKT', ...
    'currentTolerance','potentialTolerance','densityTolerance','poissonResidualTolerance','maxPotentialStepV', ...
    'chargeNeutralityTolerance'};
for k = 1:numel(positive)
    validateattributes(options.(positive{k}),{'numeric'},{'real','finite','scalar','positive'},mfilename,positive{k});
end
integers = {'maxIterations','maxEnergyRefinements','maxEnergyEvaluations','kzPoints'};
for k = 1:numel(integers)
    validateattributes(options.(integers{k}),{'numeric'},{'real','finite','scalar','integer','positive'});
end
validateattributes(options.bias,{'numeric'},{'real','finite','nonempty','vector'});
validateattributes(options.modeOrder,{'numeric'},{'real','finite','scalar','integer','nonnegative'});
validateattributes(options.andersonDepth,{'numeric'},{'real','finite','scalar','integer','nonnegative'});
validateattributes(options.mixing,{'numeric'},{'real','finite','scalar','>',0,'<=',1});
validateattributes(options.verbose,{'logical'},{'scalar'});
options.bias = options.bias(:).';
end
