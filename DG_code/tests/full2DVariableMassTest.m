classdef full2DVariableMassTest < matlab.unittest.TestCase
    %FULL2DVARIABLEMASSTEST Tests for the endpoint-mass kinetic operator.

    methods (TestClassSetup)
        function addDGPath(testCase)
            thisFile = mfilename('fullpath');
            testDir = fileparts(thisFile);
            dgDir = fileparts(testDir);
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(dgDir));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(dgDir, 'core')));
        end
    end

    methods (Test)
        function constantMassRecoversExistingTransport(testCase)
            mat = makeVariableMassMat(false);
            p = initParams_full2D(mat, mat.V, 0.2, 0.1);
            boundary = get_Boundary_full2D(mat, p, 0.2, 0.1, mat.V);

            [Aconstant, rhsConstant] = get_Diff_full2D(mat, p, boundary);
            [Avariable, rhsVariable, info] = ...
                get_Diff_variableMass_full2D(mat, p, boundary);

            testCase.verifyEqual(Avariable, Aconstant, AbsTol=1e-13);
            testCase.verifyEqual(rhsVariable, rhsConstant, AbsTol=1e-13);
            testCase.verifyEqual(info.massModel, ...
                'position-dependent-bdd');
            testCase.verifyEqual(info.mass.X.plusCorrectionRange, [0, 0], ...
                AbsTol=1e-13);
            testCase.verifyEqual(info.mass.Y.plusCorrectionRange, [0, 0], ...
                AbsTol=1e-13);
        end

        function matrixFreeMatchesAssemblyAcrossMassStep(testCase)
            setup = buildVariableMassSetup();
            u = deterministicVector(setup.p.index.nTotal);

            explicit = setup.A*u;
            matrixFree = setup.info.apply(u);
            verifyRelativeSmall(testCase, matrixFree-explicit, ...
                explicit, 2e-12);
            testCase.verifyGreaterThan(norm(setup.A-setup.referenceA, 'fro'), 0);
            testCase.verifyGreaterThan( ...
                diff(setup.info.mass.X.relativeMassRange), 0);
            testCase.verifyGreaterThan( ...
                diff(setup.info.mass.Y.relativeMassRange), 0);
            testCase.verifyLessThanOrEqual( ...
                setup.info.mass.X.endpointExchangeRelativeDefect, 1e-12);
            testCase.verifyLessThanOrEqual( ...
                setup.info.mass.Y.endpointExchangeRelativeDefect, 1e-12);
        end

        function diagonalAndRelativeBlockMatchAssembledMatrix(testCase)
            setup = buildVariableMassSetup();
            diagonal = setup.info.getDiagonal();
            expectedDiagonal = full(diag(setup.A));
            verifyRelativeSmall(testCase, diagonal-expectedDiagonal, ...
                expectedDiagonal, 2e-12);

            centerId = ceil(setup.p.dg.nCenterDof/2);
            ids = (centerId-1)*setup.p.relative.nDof ...
                + (1:setup.p.relative.nDof);
            blockData = setup.info.getRelativeBlockData();
            expectedBlock = full(setup.A(ids, ids));
            verifyRelativeSmall(testCase, ...
                full(blockData.getBlock(centerId))-expectedBlock, ...
                expectedBlock, 2e-12);
        end

        function rowAbsPreconditionerBoundCoversVariableOperator(testCase)
            setup = buildVariableMassSetup();
            rowBound = setup.info.getRowAbsSum();
            exactRowAbs = full(sum(abs(setup.A), 2));
            tolerance = 1e-12*max(1, max(exactRowAbs));

            testCase.verifyGreaterThanOrEqual(rowBound + tolerance, ...
                exactRowAbs);
            testCase.verifyTrue(all(isfinite(rowBound)));
        end

        function massCorrectionIsWeightedSkewHermitian(testCase)
            setup = buildVariableMassSetup();
            correction = setup.A-setup.referenceA;
            weightX = kron(setup.p.dg.X.jacobian(:), ...
                setup.p.dg.X.w(:));
            weightY = kron(setup.p.dg.Y.jacobian(:), ...
                setup.p.dg.Y.w(:));
            centerWeights = kron(weightY, weightX);
            weights = kron(centerWeights, ...
                ones(setup.p.relative.nDof, 1));
            W = spdiags(weights, 0, setup.p.index.nTotal, ...
                setup.p.index.nTotal);
            weightedCorrection = W*correction;
            defect = weightedCorrection + weightedCorrection';

            verifyRelativeSmall(testCase, defect, weightedCorrection, 2e-12);
        end

        function systemBuilderSelectsVariableMassPath(testCase)
            mat = makeVariableMassMat(true);
            mat.dg.params.full2D_massModel = 'position-dependent';
            mat.dg.params.full2D_assembleDriftMatrix = true;
            mat.dg.params.full2D_potentialDiscretization = 'collocated';
            p = initParams_full2D(mat, mat.V, 0.2, 0.1);

            [A, rhs, info] = get_SysM_full2D(mat, p, mat.V, 0.2, 0.1);
            u = deterministicVector(p.index.nTotal);

            testCase.verifyEqual(info.massModel, ...
                'position-dependent-bdd');
            testCase.verifyEqual(info.diff.massModel, ...
                'position-dependent-bdd');
            expected = A*u;
            verifyRelativeSmall(testCase, info.apply(u)-expected, ...
                expected, 2e-12);
            testCase.verifySize(rhs, [p.index.nTotal, 1]);
        end

        function observablesUseLocalMassOnDGGrid(testCase)
            mat = makeVariableMassMat(true);
            mat.dg.params.full2D_massModel = 'position-dependent';
            p = initParams_full2D(mat, mat.V, 0.2, 0.1);
            rho = deterministicVector(p.index.nTotal);

            [~, ~, ~, info] = get_Observables_full2D(rho, p, mat);

            testCase.verifyEqual(info.currentScale.massModel, ...
                'position-dependent-bdd');
            testCase.verifyGreaterThan( ...
                diff(info.currentScale.inverseRelativeMassXRange), 0);
            testCase.verifyGreaterThan( ...
                diff(info.currentScale.inverseRelativeMassYRange), 0);
        end
    end
end

function setup = buildVariableMassSetup()
mat = makeVariableMassMat(true);
p = initParams_full2D(mat, mat.V, 0.2, 0.1);
boundary = get_Boundary_full2D(mat, p, 0.2, 0.1, mat.V);
[referenceA, ~] = get_Diff_full2D(mat, p, boundary);
[A, rhs, info] = get_Diff_variableMass_full2D(mat, p, boundary);
setup = struct('mat', mat, 'p', p, 'boundary', boundary, ...
    'referenceA', referenceA, 'A', A, 'rhs', rhs, 'info', info);
end

function mat = makeVariableMassMat(withStep)
mat = struct;
mat.type = 'Full2DVariableMassTest';
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
mat.me_y_ch = 0.052;
mat.V = zeros(mat.Nx, mat.Ny);
mat.Nd = ones(mat.Nx, mat.Ny);

massX = mat.me_x_ch*ones(mat.Nx, mat.Ny);
massY = mat.me_y_ch*ones(mat.Nx, mat.Ny);
if withStep
    massX(3:end, 2:4) = 0.083;
    massY(2:4, [1, 5]) = 0.19;
end
mat.me_x = reshape(massX, 1, mat.Nx, mat.Ny);
mat.me_y = reshape(massY, 1, mat.Nx, mat.Ny);

mat.dg.params = struct;
mat.dg.params.full2D_coordinate_scale = 1;
mat.dg.params.N_K_X = 3;
mat.dg.params.N_K_Y = 3;
mat.dg.params.N_rho_x = 5;
mat.dg.params.N_rho_y = 5;
mat.dg.params.L_rho_x = 1;
mat.dg.params.L_rho_y = 1;
mat.dg.params.full2D_fluxType = 'rusanov';
mat.dg.params.full2D_thetaLF = 0.03;
mat.dg.params.full2D_assembleDiffMatrix = true;
mat.dg.params.full2D_reservoirModel = 'material-default';
end

function u = deterministicVector(n)
u = sin((1:n).') + 1i*cos((1:n).');
end

function verifyRelativeSmall(testCase, difference, reference, tolerance)
relativeError = norm(difference(:))/max(norm(reference(:)), eps);
testCase.verifyLessThanOrEqual(relativeError, tolerance);
end
