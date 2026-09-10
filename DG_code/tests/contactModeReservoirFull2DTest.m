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
