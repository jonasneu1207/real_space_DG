classdef contactModeReservoirFull2DTest < matlab.unittest.TestCase
    %CONTACTMODERESERVOIRFULL2DTEST Tests contact-mode Source/Drain data.

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
        function contactModeReservoirHasExpectedShapeAndOrigin(testCase)
            mat = makeContactModeMat();
            p = initParams_full2D(mat, mat.V, 2.0, 1.8);

            [inflow, info] = get_InflowBoundary_full2D(mat, p, 2.0, 1.8, mat.V);

            testCase.verifyEqual(info.source.origin, 'contact-modes');
            testCase.verifyEqual(info.source.reservoirModel, 'contact-modes');
            testCase.verifySize(inflow.source.rhoBoundary, ...
                [p.relative.nDof, p.dg.Y.nDof]);
            testCase.verifyGreaterThan(norm(inflow.source.rhoBoundary, 'fro'), 0);
            testCase.verifyEqual(info.source.defaultReservoir.modeCount, 2);
            testCase.verifyGreaterThan(info.source.defaultReservoir.contact.HNnz, 0);
        end

        function contactModeReservoirVariesAlongY(testCase)
            mat = makeContactModeMat();
            p = initParams_full2D(mat, mat.V, 2.0, 1.8);

            [inflow, ~] = get_InflowBoundary_full2D(mat, p, 2.0, 1.8, mat.V);
            columnNorms = vecnorm(inflow.source.rhoBoundary, 2, 1);

            testCase.verifyGreaterThan(max(columnNorms) - min(columnNorms), 0);
        end

        function contactModesUseContinuumYNormalization(testCase)
            mat = makeContactModeMat();
            p = initParams_full2D(mat, mat.V, 2.0, 1.8);

            [~, info] = get_InflowBoundary_full2D(mat, p, 2.0, 1.8, mat.V);
            normalizedIntegral = ...
                info.source.defaultReservoir.modeNormalization.normalizedIntegral;

            testCase.verifyEqual(normalizedIntegral, ...
                ones(size(normalizedIntegral)), AbsTol=1e-10);
        end

        function contactModeKzCutoffCoversOccupiedModesAndThermalTail(testCase)
            mat = makeContactModeMat();
            mat.dg.params.full2D_contactModeIntegrateKz = true;
            mat.dg.params.full2D_contactMode_Nkz = 31;
            mat.dg.params.full2D_contactModeKzEnergyWindow = 0.05;
            mat.dg.params.full2D_contactModeKzThermalTail = 8;
            Ef = 2.0;
            p = initParams_full2D(mat, mat.V, Ef, Ef);

            [~, info] = get_InflowBoundary_full2D(mat, p, Ef, Ef, mat.V);
            reservoir = info.source.defaultReservoir;
            kz = reservoir.kz;
            thermalTailEV = 8*1.38064852e-23*mat.Temp/1.602176634e-19;
            expectedAdaptiveWindow = max(Ef-min(reservoir.modeEnergy), 0) ...
                + thermalTailEV;

            testCase.verifyTrue(kz.integrate);
            testCase.verifyEqual(kz.minimumEnergyWindowEV, 0.05);
            testCase.verifyEqual(kz.adaptiveEnergyWindowEV, ...
                expectedAdaptiveWindow, RelTol=1e-12);
            testCase.verifyEqual(kz.energyWindowEV, ...
                max(0.05, expectedAdaptiveWindow), RelTol=1e-12);
            testCase.verifyEqual(kz.thermalTailKBT, 8);
            testCase.verifyEqual(kz.nKz, 31);
        end

        function contactModeFermiLevelIsSolvedFromContactNeutrality(testCase)
            mat = makeContactModeMat();
            mat.N_s = 2e18;
            mat.N_d = 2e18;
            mat.W_c = 4;
            mat.dg.params.full2D_contactModeIntegrateKz = true;
            mat.dg.params.full2D_contactMode_Nkz = 101;
            mat.dg.params.full2D_contactModeFermiModel = 'contact-neutrality';
            inputEfL = 5.0;
            inputEfR = -2.0;
            p = initParams_full2D(mat, mat.V, inputEfL, inputEfR);

            [inflow, info] = get_InflowBoundary_full2D( ...
                mat, p, inputEfL, inputEfR, mat.V);
            source = info.source.defaultReservoir;
            drain = info.drain.defaultReservoir;
            expectedTarget = mat.N_s*1e6*mat.W_c*1e-9;

            testCase.verifyEqual(source.fermiLevel.model, 'contact-neutrality');
            testCase.verifyEqual(source.fermiLevel.targetSheetDensityM2, ...
                expectedTarget, RelTol=1e-14);
            testCase.verifyLessThan(source.fermiLevel.relativeResidual, 1e-9);
            testCase.verifyLessThan(drain.fermiLevel.relativeResidual, 1e-9);
            testCase.verifyNotEqual(source.fermiLevelUsed, inputEfL);
            testCase.verifyNotEqual(drain.fermiLevelUsed, inputEfR);
            testCase.verifyEqual(source.fermiLevelUsed, ...
                drain.fermiLevelUsed, AbsTol=1e-10);
            testCase.verifyEqual(inflow.source.fermiLevel, ...
                source.fermiLevelUsed, AbsTol=1e-14);
            testCase.verifyEqual(inflow.drain.fermiLevel, ...
                drain.fermiLevelUsed, AbsTol=1e-14);
        end

        function externalContactModeFermiLevelRemainsSelectable(testCase)
            mat = makeContactModeMat();
            mat.N_s = 2e18;
            mat.N_d = 2e18;
            mat.W_c = 4;
            mat.dg.params.full2D_contactModeIntegrateKz = true;
            mat.dg.params.full2D_contactMode_Nkz = 31;
            mat.dg.params.full2D_contactModeFermiModel = 'external';
            Ef = 1.25;
            p = initParams_full2D(mat, mat.V, Ef, Ef);

            [inflow, info] = get_InflowBoundary_full2D(mat, p, Ef, Ef, mat.V);

            testCase.verifyEqual( ...
                info.source.defaultReservoir.fermiLevel.model, 'external');
            testCase.verifyEqual(inflow.source.fermiLevel, Ef);
            testCase.verifyEqual( ...
                info.source.defaultReservoir.fermiLevelUsed, Ef);
        end

        function equilibriumContactFermiHelperUsesSameModeModel(testCase)
            mat = makeContactModeMat();
            mat.N_s = 2e18;
            mat.N_d = 2e18;
            mat.W_c = 4;
            mat.dg.params.full2D_contactModeIntegrateKz = true;
            mat.dg.params.full2D_contactMode_Nkz = 31;
            mat.dg.params.full2D_contactModeFermiModel = 'contact-neutrality';

            [EfL, EfR, info] = solve_fermi_contactModes_full2D( ...
                mat, mat.V, 5.0, -2.0);

            testCase.verifyEqual(EfL, info.source.usedEV, AbsTol=1e-14);
            testCase.verifyEqual(EfR, info.drain.usedEV, AbsTol=1e-14);
            testCase.verifyEqual(EfL, EfR, AbsTol=1e-10);
            testCase.verifyTrue(info.source.converged);
            testCase.verifyTrue(info.drain.converged);
        end

        function contactModeReservoirDiffersFromMaterialDefault(testCase)
            contactMode = makeContactModeMat();
            materialDefault = contactMode;
            materialDefault.dg.params.full2D_reservoirModel = 'material-default';
            pContact = initParams_full2D(contactMode, contactMode.V, 2.0, 1.8);
            pMaterial = initParams_full2D(materialDefault, materialDefault.V, 2.0, 1.8);

            [contactInflow, ~] = get_InflowBoundary_full2D( ...
                contactMode, pContact, 2.0, 1.8, contactMode.V);
            [materialInflow, ~] = get_InflowBoundary_full2D( ...
                materialDefault, pMaterial, 2.0, 1.8, materialDefault.V);

            difference = norm(contactInflow.source.rhoBoundary ...
                - materialInflow.source.rhoBoundary, 'fro');

            testCase.verifyGreaterThan(difference, 0);
        end

        function explicitRhoDataOverridesContactModeReservoir(testCase)
            mat = makeContactModeMat();
            p = initParams_full2D(mat, mat.V, 2.0, 1.8);
            expectedSource = ones(p.relative.nDof, p.dg.Y.nDof);
            expectedDrain = 2*ones(p.relative.nDof, p.dg.Y.nDof);
            mat.dg.params.full2D_sourceRho = expectedSource;
            mat.dg.params.full2D_drainRho = expectedDrain;

            [inflow, info] = get_InflowBoundary_full2D(mat, p, 2.0, 1.8, mat.V);

            testCase.verifyEqual(info.source.origin, 'full2D_sourceRho');
            testCase.verifyEqual(info.drain.origin, 'full2D_drainRho');
            testCase.verifyEqual(inflow.source.rhoBoundary, expectedSource);
            testCase.verifyEqual(inflow.drain.rhoBoundary, expectedDrain);
        end

        function contactModeReservoirRespectsContactMask(testCase)
            mat = makeMaskedContactModeMat();
            p = initParams_full2D(mat, mat.V, 2.0, 1.8);
            expectedMask = get_ContactMask_full2D( ...
                mat, p, 'source', p.dg.Y.nodes(:));

            [inflow, ~] = get_InflowBoundary_full2D(mat, p, 2.0, 1.8, mat.V);

            testCase.verifyEqual(inflow.source.contactMask, expectedMask);
            testCase.verifyGreaterThan( ...
                norm(inflow.source.rhoBoundary(:, expectedMask), 'fro'), 0);
            testCase.verifyEqual(inflow.source.rhoBoundary(:, ~expectedMask), ...
                zeros(p.relative.nDof, nnz(~expectedMask)), AbsTol=1e-12);
        end
    end
end

function mat = makeContactModeMat()
mat = struct;
mat.type = 'ContactModeReservoirFull2DTest';
mat.Nx = 5;
mat.Ny = 9;
mat.x = linspace(0, 4, mat.Nx);
mat.y = linspace(-2, 2, mat.Ny);
mat.dx = mat.x(2) - mat.x(1);
mat.dy = mat.y(2) - mat.y(1);
mat.Temp = 300;
mat.deg_factor = 1;
mat.Eg_c = 0.74;
mat.Eg_ox = 8.8;
mat.me_x_ch = 1.0;
mat.me_y_ch = 1.0;
mat.me_z_ch = 1.0;
mat.n_of_modes = 2;

[~, Y] = ndgrid(mat.x, mat.y);
mat.V = 0.05 + 0.08*Y.^2;
mat.V(:, [1, end]) = 1.0;
mat.Nd = ones(mat.Nx, mat.Ny);
mat.me_x = ones(1, mat.Nx, mat.Ny);
mat.me_y = ones(1, mat.Nx, mat.Ny);
mat.me_z = ones(1, mat.Nx, mat.Ny);

mat.dg.params = struct;
mat.dg.params.full2D_coordinate_scale = 1e-9;
mat.dg.params.N_K_chi = 3;
mat.dg.params.N_K_X = 3;
mat.dg.params.N_K_Y = 3;
mat.dg.params.N_rho_x = 5;
mat.dg.params.N_rho_y = 5;
mat.dg.params.L_rho_x = 8e-9;
mat.dg.params.L_rho_y = 4e-9;
mat.dg.params.full2D_reservoirModel = 'contact-modes';
mat.dg.params.full2D_contactModeCount = 2;
mat.dg.params.full2D_contactModeIntegrateKz = false;
mat.dg.params.full2D_contactMode_Nkx = 7;
end

function mat = makeMaskedContactModeMat()
mat = makeContactModeMat();
mat.Nd = zeros(mat.Nx, mat.Ny);
contactY = mat.y >= -1 & mat.y <= 1;
mat.Nd(1, contactY) = 1;
mat.Nd(end, contactY) = 1;
end
