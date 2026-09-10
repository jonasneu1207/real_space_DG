classdef full2DDiffOperatorTest < matlab.unittest.TestCase
    %FULL2DDIFFOPERATORTEST Smoke tests for the Full-2D rho transport matrix.

    methods (TestClassSetup)
        function addDGPath(testCase)
            thisFile = mfilename('fullpath');
            testDir = fileparts(thisFile);
            dgDir = fileparts(testDir);
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(dgDir));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(dgDir, 'core')));
        end
    end

    methods (Test)
        function rusanovAssemblyIsSparseAndSized(testCase)
            setup = buildDiffSetup('rusanov', true);

            testCase.verifyTrue(issparse(setup.A));
            testCase.verifySize(setup.A, [setup.p.index.nTotal, setup.p.index.nTotal]);
            testCase.verifySize(setup.rhs, [setup.p.index.nTotal, 1]);
            testCase.verifyTrue(setup.info.fullMatrixAssembled);
            testCase.verifyEqual(setup.info.fluxType, 'rusanov');
            testCase.verifyGreaterThan(nnz(setup.A), 0);
            testCase.verifyGreaterThan(norm(setup.rhs), 0);
            testCase.verifyGreaterThan(setup.info.betaLF(1), 0);
            testCase.verifyGreaterThan(setup.info.betaLF(2), 0);
        end

        function rusanovAddsInteriorStabilizationToCentral(testCase)
            central = buildDiffSetup('central', true);
            rusanov = buildDiffSetup('rusanov', true);

            testCase.verifyEqual(central.info.betaLF, [0, 0], AbsTol=1e-14);
            testCase.verifyGreaterThan(rusanov.info.betaLF(1), 0);
            testCase.verifyGreaterThan(rusanov.info.betaLF(2), 0);
            testCase.verifyGreaterThan(norm(rusanov.A - central.A, 'fro'), 0);
        end

        function rhoFluxParameterIsAcceptedAsFallback(testCase)
            setup = buildRhoFluxFallbackSetup();

            testCase.verifyEqual(setup.info.fluxType, 'central');
            testCase.verifyEqual(setup.info.betaLF, [0, 0], AbsTol=1e-14);
        end

        function matrixFreeApplyMatchesAssembledMatrix(testCase)
            setup = buildDiffSetup('rusanov', true);
            u = deterministicVector(setup.p.index.nTotal);

            explicitResult = setup.A*u;
            matrixFreeResult = setup.info.apply(u);

            testCase.verifyEqual(matrixFreeResult, explicitResult, AbsTol=1e-10);
        end

        function autoModeCanKeepMatrixFree(testCase)
            setup = buildMatrixFreeSetup();
            u = deterministicVector(setup.p.index.nTotal);

            testCase.verifyEmpty(setup.A);
            testCase.verifyEmpty(setup.rhs);
            testCase.verifyFalse(setup.info.fullMatrixAssembled);
            testCase.verifyTrue(setup.info.matrixFreeAvailable);
            testCase.verifySize(setup.info.apply(u), [setup.p.index.nTotal, 1]);
        end

        function sourceRhsUsesOnlyContactMask(testCase)
            setup = buildMaskedSourceSetup();
            rhsByCenterColumn = reshape(setup.rhs, ...
                setup.p.relative.nDof, setup.p.dg.nCenterDof);

            testCase.verifyEqual(full(rhsByCenterColumn(:, setup.closedRows)), ...
                zeros(setup.p.relative.nDof, numel(setup.closedRows)), AbsTol=1e-12);
            testCase.verifyGreaterThan( ...
                norm(full(rhsByCenterColumn(:, setup.contactRows)), 'fro'), 0);
        end
    end
end

function setup = buildRhoFluxFallbackSetup()
mat = makeSmallDiffMat();
mat.dg.params.rho_flux = 'central';
mat.dg.params.full2D_assembleDiffMatrix = true;
if isfield(mat.dg.params, 'full2D_fluxType')
    mat.dg.params = rmfield(mat.dg.params, 'full2D_fluxType');
end
p = initParams_full2D(mat, mat.V, 0.2, 0.1);
boundary = get_Boundary_full2D(mat, p, 0.2, 0.1, mat.V);
[A, rhs, info] = get_Diff_full2D(mat, p, boundary);

setup = struct;
setup.p = p;
setup.A = A;
setup.rhs = rhs;
setup.info = info;
end

function setup = buildDiffSetup(fluxType, assembleMatrix)
mat = makeSmallDiffMat();
mat.dg.params.full2D_fluxType = fluxType;
mat.dg.params.full2D_assembleDiffMatrix = assembleMatrix;
p = initParams_full2D(mat, mat.V, 0.2, 0.1);
boundary = get_Boundary_full2D(mat, p, 0.2, 0.1, mat.V);
[A, rhs, info] = get_Diff_full2D(mat, p, boundary);

setup = struct;
setup.mat = mat;
setup.p = p;
setup.boundary = boundary;
setup.A = A;
setup.rhs = rhs;
setup.info = info;
end

function setup = buildMatrixFreeSetup()
mat = makeSmallDiffMat();
mat.dg.params.full2D_fluxType = 'rusanov';
mat.dg.params.full2D_assembleDiffMatrix = 'auto';
mat.dg.params.full2D_maxAssembledDof = 100;
p = initParams_full2D(mat, mat.V, 0.2, 0.1);
boundary = get_Boundary_full2D(mat, p, 0.2, 0.1, mat.V);
[A, rhs, info] = get_Diff_full2D(mat, p, boundary);

setup = struct;
setup.p = p;
setup.A = A;
setup.rhs = rhs;
setup.info = info;
end

function setup = buildMaskedSourceSetup()
mat = makeSmallDiffMat();
p = initParams_full2D(mat, mat.V, 0.2, 0.1);
sourceProfile = repmat((1:p.relative.nDof).', 1, p.dg.Y.nDof);
contactMask = false(p.dg.Y.nDof, 1);
contactMask(2:4) = true;
mat.dg.params.full2D_fluxType = 'rusanov';
mat.dg.params.full2D_assembleDiffMatrix = true;
mat.dg.params.full2D_sourceRho = sourceProfile;
mat.dg.params.full2D_drainRho = zeros(p.relative.nDof, p.dg.Y.nDof);
mat.dg.params.full2D_sourceContactMask = contactMask;
mat.dg.params.full2D_drainContactMask = false(p.dg.Y.nDof, 1);

boundary = get_Boundary_full2D(mat, p, 0.2, 0.1, mat.V);
[~, rhs, ~] = get_Diff_full2D(mat, p, boundary);

setup = struct;
setup.p = p;
setup.rhs = rhs;
setup.closedRows = centerRowsForYMask(p, ~contactMask);
setup.contactRows = centerRowsForYMask(p, contactMask);
end

function rows = centerRowsForYMask(p, yMask)
xRows = (1:p.dg.X.nLocal).';
yIds = find(yMask(:)).';
rows = zeros(numel(xRows)*numel(yIds), 1);
cursor = 1;
for iy = yIds
    ids = sub2ind([p.dg.X.nDof, p.dg.Y.nDof], xRows, iy*ones(size(xRows)));
    rows(cursor:cursor+numel(ids)-1) = ids;
    cursor = cursor + numel(ids);
end
end

function mat = makeSmallDiffMat()
mat = struct;
mat.type = 'Full2DDiffTest';
mat.Nx = 5;
mat.Ny = 5;
mat.x = linspace(0, 1, mat.Nx);
mat.y = linspace(-0.5, 0.5, mat.Ny);
mat.dx = mat.x(2)-mat.x(1);
mat.dy = mat.y(2)-mat.y(1);
mat.Temp = 300;
mat.deg_factor = 1;
mat.Eg_c = 0.74;
mat.Eg_ox = 8.8;
mat.me_x_ch = 0.041;
mat.me_y_ch = 0.041;
mat.V = zeros(mat.Nx, mat.Ny);
mat.Nd = ones(mat.Nx, mat.Ny);
mat.me_x = 0.041*ones(1, mat.Nx, mat.Ny);
mat.me_y = 0.041*ones(1, mat.Nx, mat.Ny);
mat.dg.params = struct;
mat.dg.params.full2D_coordinate_scale = 1;
mat.dg.params.N_K_chi = 3;
mat.dg.params.N_K_X = 3;
mat.dg.params.N_K_Y = 3;
mat.dg.params.N_rho_x = 5;
mat.dg.params.N_rho_y = 5;
mat.dg.params.L_rho_x = 1;
mat.dg.params.L_rho_y = 1;
mat.dg.params.full2D_thetaLF = 0.03;
end

function u = deterministicVector(n)
u = reshape(sin((1:n).') + 1i*cos((1:n).'), [], 1);
end
