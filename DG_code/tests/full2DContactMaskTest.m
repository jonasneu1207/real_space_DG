classdef full2DContactMaskTest < matlab.unittest.TestCase
    %FULL2DCONTACTMASKTEST Tests masked Source/Drain openings on X faces.

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
        function sourceMaskUsesDopedBuildDeviceEdge(testCase)
            mat = makeMaskedMat();
            p = initParams_full2D(mat, mat.V, 0.2, 0.1);

            [mask, info] = get_ContactMask_full2D(mat, p, 'source', p.dg.Y.nodes(:));
            expected = p.dg.Y.nodes(:) >= 0 & p.dg.Y.nodes(:) <= 1;

            testCase.verifyEqual(mask, expected);
            testCase.verifyEqual(info.origin, 'buildDevice-Nd-edge');
            testCase.verifyEqual(info.nClosedDof, nnz(~expected));
        end

        function sourceBoundaryZerosReservoirOutsideMask(testCase)
            mat = makeMaskedMat();
            p = initParams_full2D(mat, mat.V, 0.2, 0.1);
            mat.dg.params.full2D_sourceRho = ones(p.relative.nDof, p.dg.Y.nDof);
            mat.dg.params.full2D_drainRho = 2*ones(p.relative.nDof, p.dg.Y.nDof);

            [inflow, ~] = get_InflowBoundary_full2D(mat, p, 0.2, 0.1, mat.V);
            mask = p.dg.Y.nodes(:) >= 0 & p.dg.Y.nodes(:) <= 1;

            testCase.verifyEqual(inflow.source.contactMask, mask);
            testCase.verifyEqual(inflow.source.rhoBoundary(:, mask), ...
                ones(p.relative.nDof, nnz(mask)));
            testCase.verifyEqual(inflow.source.rhoBoundary(:, ~mask), ...
                zeros(p.relative.nDof, nnz(~mask)));
        end

        function maskedFluxUsesSpecularXReflection(testCase)
            [mat, p] = makeToyMaskedFluxProblem();
            mat.dg.params.full2D_sourceRho = repmat([10; 20; 30; 40], 1, 2);
            mat.dg.params.full2D_drainRho = repmat([100; 200; 300; 400], 1, 2);
            mat.dg.params.full2D_sourceContactMask = [true; false];
            mat.dg.params.full2D_drainContactMask = [true; false];

            [inflow, ~] = get_InflowBoundary_full2D(mat, p, 0.1, 0.0);
            rhoInside = repmat([1; 2; 3; 4], 1, 2);
            flux = apply_InflowFlux_full2D(inflow.source, rhoInside);
            expected = [[2; 2; -30; -80], [2; 2; -2; -2]];

            testCase.verifyEqual(flux, expected, AbsTol=1e-12);
            testCase.verifyEqual(inflow.source.nContactFaceDof, 1);
            testCase.verifyEqual(inflow.source.nClosedFaceDof, 1);
        end
    end
end

function mat = makeMaskedMat()
mat = struct;
mat.type = 'Full2DContactMaskTest';
mat.Nx = 5;
mat.Ny = 7;
mat.x = linspace(0, 1, mat.Nx);
mat.y = linspace(-1, 2, mat.Ny);
mat.dx = mat.x(2)-mat.x(1);
mat.dy = mat.y(2)-mat.y(1);
mat.Temp = 300;
mat.deg_factor = 1;
mat.Eg_c = 0.74;
mat.Eg_ox = 8.8;
mat.me_x_ch = 0.041;
mat.me_y_ch = 0.041;
mat.V = zeros(mat.Nx, mat.Ny);
mat.Nd = zeros(mat.Nx, mat.Ny);
contactY = mat.y >= 0 & mat.y <= 1;
mat.Nd(1, contactY) = 2e19;
mat.Nd(end, contactY) = 2e19;
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
end

function [mat, p] = makeToyMaskedFluxProblem()
mat = struct;
mat.dg.params = struct;
p = struct;
p.domain.normals.XLeft = [-1, 0];
p.domain.normals.XRight = [1, 0];
p.domain.normals.YBottom = [0, -1];
p.domain.normals.YTop = [0, 1];
p.domain.normalConvention = ...
    'Normals point out of the physical X-Y domain; negative eigenvalues of A_n are inflow.';
p.index.order = {'rho_x', 'rho_y', 'X-DG-DOF', 'Y-DG-DOF'};
p.relative.NrhoX = 4;
p.relative.NrhoY = 1;
p.relative.size = [4, 1];
p.relative.nDof = 4;
p.relative.Ax1D = spdiags([-2; -1; 1; 2], 0, 4, 4);
p.relative.Ay1D = sparse(1, 1);
p.relative.Ax = p.relative.Ax1D;
p.relative.Ay = sparse(4, 4);
p.relative.identity = speye(4);
p.dg.Y.nDof = 2;
p.dg.Y.nodes = [0; 1];
p.dg.faces.global.left.centerDofs = (1:2).';
p.dg.faces.global.right.centerDofs = (1:2).';
end
