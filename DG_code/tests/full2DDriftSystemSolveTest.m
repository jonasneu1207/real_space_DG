classdef full2DDriftSystemSolveTest < matlab.unittest.TestCase
    %FULL2DDRIFTSYSTEMSOLVETEST Smoke tests for Full-2D drift/system/solve.

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
        function driftIsSparseDiagonalAndMatrixFreeMatches(testCase)
            mat = makeSmallMat();
            mat.dg.params.full2D_assembleDriftMatrix = true;
            p = initParams_full2D(mat, mat.V, 0.2, 0.1);

            [G, info] = get_Drift_full2D(mat, p, mat.V);
            u = deterministicVector(p.index.nTotal);

            testCase.verifyTrue(issparse(G));
            testCase.verifyEqual(nnz(G - spdiags(diag(G), 0, ...
                p.index.nTotal, p.index.nTotal)), 0);
            verifyRelativeSmall(testCase, info.apply(u) - G*u, G*u, 1e-13);
            testCase.verifyTrue(info.fullPotentialAssembled);
            testCase.verifyTrue(info.matrixFreeAvailable);
        end

        function zeroRelativeSeparationHasNoPotentialDifference(testCase)
            mat = makeSmallMat();
            mat.dg.params.full2D_assembleDriftMatrix = true;
            p = initParams_full2D(mat, mat.V, 0.2, 0.1);

            [~, info] = get_Drift_full2D(mat, p, mat.V);
            diagonal = info.getDiagonal();
            [~, irx0] = min(abs(p.relative.rhoX.cells));
            [~, iry0] = min(abs(p.relative.rhoY.cells));
            relativeId = irx0 + (iry0-1)*p.relative.NrhoX;
            ids = relativeId + (0:p.dg.nCenterDof-1)*p.relative.nDof;

            testCase.verifyEqual(diagonal(ids), ...
                zeros(p.dg.nCenterDof, 1), AbsTol=1e-10);
        end

        function legacyRhoDriftScaleDoesNotAffectFull2D(testCase)
            reference = makeSmallMat();
            reference.dg.params = rmfield(reference.dg.params, ...
                'full2D_drift_scale');
            pReference = initParams_full2D(reference, reference.V, 0.2, 0.1);
            [~, referenceInfo] = get_Drift_full2D(reference, ...
                pReference, reference.V);

            legacy = reference;
            legacy.dg.params.rho_drift_scale = 1e9;
            pLegacy = initParams_full2D(legacy, legacy.V, 0.2, 0.1);
            [~, legacyInfo] = get_Drift_full2D(legacy, pLegacy, legacy.V);

            explicit = reference;
            explicit.dg.params.full2D_drift_scale = 2;
            pExplicit = initParams_full2D(explicit, explicit.V, 0.2, 0.1);
            [~, explicitInfo] = get_Drift_full2D(explicit, ...
                pExplicit, explicit.V);

            testCase.verifyEqual(legacyInfo.getDiagonal(), ...
                referenceInfo.getDiagonal(), AbsTol=1e-12);
            testCase.verifyEqual(explicitInfo.getDiagonal(), ...
                2*referenceInfo.getDiagonal(), AbsTol=1e-12);
            testCase.verifyTrue(isfield(legacyInfo.scale, ...
                'legacyRhoDriftScaleIgnored'));
        end

        function matrixFreeDriftCanSkipStoredDiagonal(testCase)
            mat = makeSmallMat();
            mat.dg.params.full2D_assembleDriftMatrix = false;
            mat.dg.params.full2D_storeDriftDiagonal = false;
            p = initParams_full2D(mat, mat.V, 0.2, 0.1);

            [G, info] = get_Drift_full2D(mat, p, mat.V);
            u = deterministicVector(p.index.nTotal);

            testCase.verifyEmpty(G);
            testCase.verifyFalse(info.fullPotentialAssembled);
            testCase.verifyFalse(info.diagonalStored);
            testCase.verifySize(info.apply(u), [p.index.nTotal, 1]);
        end

        function systemMatrixCombinesDiffAndDrift(testCase)
            mat = makeSmallMat();
            mat.dg.params.full2D_assembleDiffMatrix = true;
            mat.dg.params.full2D_assembleDriftMatrix = true;
            p = initParams_full2D(mat, mat.V, 0.2, 0.1);

            [A, rhs, info] = get_SysM_full2D(mat, p, mat.V, 0.2, 0.1);
            u = deterministicVector(p.index.nTotal);
            expected = info.diff.assemble() + info.drift.assemble();

            testCase.verifyTrue(info.fullMatrixAssembled);
            testCase.verifySize(A, [p.index.nTotal, p.index.nTotal]);
            testCase.verifySize(rhs, [p.index.nTotal, 1]);
            verifyRelativeSmall(testCase, A*u - expected*u, expected*u, 1e-13);
            verifyRelativeSmall(testCase, info.apply(u) - A*u, A*u, 1e-13);
        end

        function systemDiagonalMatchesAssembledMatrix(testCase)
            mat = makeSmallMat();
            mat.dg.params.full2D_assembleDiffMatrix = true;
            mat.dg.params.full2D_assembleDriftMatrix = true;
            p = initParams_full2D(mat, mat.V, 0.2, 0.1);

            [A, ~, info] = get_SysM_full2D(mat, p, mat.V, 0.2, 0.1);
            diagonal = info.getDiagonal();

            testCase.verifyEqual(diagonal, full(diag(A)), AbsTol=1e-12);
        end

        function autoLargeSystemKeepsSystemMatrixFree(testCase)
            mat = makeSmallMat();
            mat.dg.params.full2D_assembleDiffMatrix = 'auto';
            mat.dg.params.full2D_assembleDriftMatrix = 'auto';
            mat.dg.params.full2D_maxAssembledDof = 100;
            p = initParams_full2D(mat, mat.V, 0.2, 0.1);

            [A, rhs, info] = get_SysM_full2D(mat, p, mat.V, 0.2, 0.1);
            u = deterministicVector(p.index.nTotal);

            testCase.verifyEmpty(A);
            testCase.verifyEmpty(rhs);
            testCase.verifyFalse(info.fullMatrixAssembled);
            testCase.verifyTrue(info.matrixFreeAvailable);
            testCase.verifySize(info.apply(u), [p.index.nTotal, 1]);
        end

        function smallSolverReturnsDensityAndCurrentShapes(testCase)
            mat = makeSmallMat();
            mat.dg.params.full2D_assembleDiffMatrix = true;
            mat.dg.params.full2D_assembleDriftMatrix = true;
            mat.dg.params.full2D_solve = true;
            mat.dg.params.full2D_solver = 'direct';

            DG = solve_transport_DG_full2D(mat, mat.V, 0.2, 0.1);

            testCase.verifyEqual(DG.status, 'solved-direct');
            testCase.verifySize(DG.rho, [DG.p.index.nTotal, 1]);
            testCase.verifySize(DG.n, [DG.p.dg.X.nDof, mat.Ny]);
            testCase.verifySize(DG.jx, [DG.p.dg.X.nDof, mat.Ny]);
            testCase.verifySize(DG.j, [DG.p.dg.X.nDof, 1]);
            residual = norm(DG.info.apply(DG.rho) - DG.rhs)/max(norm(DG.rhs), eps);
            testCase.verifyLessThan(residual, 1e-9);

            [slice, axes, sliceInfo] = DG.getRhoSlice2D({'X', 'rho_x'}, ...
                struct('YIndex', 2, 'rho_yIndex', 3));
            rho4D = reshape(DG.rho, DG.p.index.arraySize);
            expected = squeeze(rho4D(:, 3, :, 2)).';

            testCase.verifySize(slice, [DG.p.dg.X.nDof, DG.p.relative.NrhoX]);
            testCase.verifyEqual(slice, expected, AbsTol=1e-12);
            testCase.verifyEqual(axes.freeDims, {'X', 'rho_x'});
            testCase.verifyEqual(sliceInfo.fixedDims, {'rho_y', 'Y'});
            testCase.verifyEqual(sliceInfo.indices.Y.index, 2);
            testCase.verifyEqual(sliceInfo.indices.rho_y.index, 3);
        end

        function rhoSlicePreservesRequestedDimensionOrder(testCase)
            mat = makeSmallMat();
            p = initParams_full2D(mat, mat.V, 0.2, 0.1);
            rho = (1:p.index.nTotal).';

            [slice, axes, sliceInfo] = get_RhoSlice2D_full2D(rho, p, ...
                {'rho_y', 'X'}, struct('YIndex', 2, 'rho_xIndex', 3));
            rho4D = reshape(rho, p.index.arraySize);
            expected = squeeze(rho4D(3, :, :, 2));

            testCase.verifySize(slice, [p.relative.NrhoY, p.dg.X.nDof]);
            testCase.verifyEqual(slice, expected);
            testCase.verifyEqual(axes.freeDims, {'rho_y', 'X'});
            testCase.verifyEqual(sliceInfo.indices.Y.index, 2);
            testCase.verifyEqual(sliceInfo.indices.rho_x.index, 3);
        end

        function bicgstabCanUseJacobiPreconditioner(testCase)
            mat = makeSmallMat();
            directMat = mat;
            directMat.dg.params.full2D_assembleDiffMatrix = true;
            directMat.dg.params.full2D_assembleDriftMatrix = true;
            directMat.dg.params.full2D_solve = true;
            directMat.dg.params.full2D_solver = 'direct';

            reference = solve_transport_DG_full2D(directMat, ...
                directMat.V, 0.2, 0.1);

            mat.dg.params.full2D_assembleDiffMatrix = false;
            mat.dg.params.full2D_assembleDriftMatrix = false;
            mat.dg.params.full2D_solve = true;
            mat.dg.params.full2D_solver = 'bicgstab';
            mat.dg.params.full2D_preconditioner = 'jacobi';
            mat.dg.params.full2D_solverTol = 1e-10;
            mat.dg.params.full2D_solverMaxIt = 2000;

            DG = solve_transport_DG_full2D(mat, mat.V, 0.2, 0.1);

            testCase.verifyTrue(DG.info.solve.preconditioner.enabled);
            testCase.verifyEqual(DG.info.solve.preconditioner.requested, ...
                'jacobi');
            testCase.verifyEqual(DG.info.solve.preconditioner.method, ...
                'rowabs');
            testCase.verifyTrue(isfield(DG.info.solve.preconditioner, ...
                'jacobiDiagonal'));
            testCase.verifyNotEmpty(DG.rho);
            testCase.verifyEqual(DG.info.solve.flag, 0);
            testCase.verifyTrue(isfinite(DG.info.solve.relres));
            verifyRelativeSmall(testCase, DG.rho - reference.rho, ...
                reference.rho, 1e-7);
        end

        function gmresCanUseRhoBlockJacobiPreconditioner(testCase)
            mat = makeSmallMat();
            mat.dg.params.N_rho_x = 4;
            mat.dg.params.N_rho_y = 4;

            mat.dg.params.full2D_assembleDiffMatrix = false;
            mat.dg.params.full2D_assembleDriftMatrix = false;
            mat.dg.params.full2D_solve = true;
            mat.dg.params.full2D_solver = 'gmres';
            mat.dg.params.full2D_preconditioner = 'blockjacobi';
            mat.dg.params.full2D_solverTol = 1e-10;
            mat.dg.params.full2D_solverMaxIt = 500;

            DG = solve_transport_DG_full2D(mat, mat.V, 0.2, 0.1);
            residual = norm(DG.info.apply(DG.rho) - DG.rhs) ...
                /max(norm(DG.rhs), eps);

            testCase.verifyTrue(DG.info.solve.preconditioner.enabled);
            testCase.verifyEqual(DG.info.solve.preconditioner.requested, ...
                'blockjacobi');
            testCase.verifyEqual(DG.info.solve.preconditioner.method, ...
                'blockjacobi-rho');
            testCase.verifyEqual(DG.info.solve.preconditioner.blockDof, ...
                DG.p.relative.nDof);
            testCase.verifyEqual(DG.info.solve.preconditioner.nBlocks, ...
                DG.p.dg.nCenterDof);
            testCase.verifyNotEmpty(DG.rho);
            testCase.verifyEqual(DG.info.solve.flag, 0);
            testCase.verifyTrue(isfinite(DG.info.solve.relres));
            testCase.verifyLessThan(residual, 1e-7);
        end

        function blockJacobiFallsBackForSingularRhoBlocks(testCase)
            mat = makeSmallMat();
            mat.dg.params.full2D_assembleDiffMatrix = false;
            mat.dg.params.full2D_assembleDriftMatrix = false;
            mat.dg.params.full2D_solve = true;
            mat.dg.params.full2D_solver = 'bicgstab';
            mat.dg.params.full2D_preconditioner = 'blockjacobi';
            mat.dg.params.full2D_solverTol = 1e-10;
            mat.dg.params.full2D_solverMaxIt = 2000;

            DG = solve_transport_DG_full2D(mat, mat.V, 0.2, 0.1);

            testCase.verifyTrue(DG.info.solve.preconditioner.enabled);
            testCase.verifyEqual(DG.info.solve.preconditioner.requested, ...
                'blockjacobi');
            testCase.verifyEqual(DG.info.solve.preconditioner.method, ...
                'rowabs');
            testCase.verifyTrue(isfield(DG.info.solve.preconditioner, ...
                'blockJacobiError'));
            testCase.verifyEqual(DG.info.solve.flag, 0);
        end

        function bicgstabRetriesRowAbsAfterBlockJacobiFailure(testCase)
            mat = makeSmallMat();
            mat.dg.params.N_rho_x = 4;
            mat.dg.params.N_rho_y = 4;
            mat.dg.params.full2D_assembleDiffMatrix = false;
            mat.dg.params.full2D_assembleDriftMatrix = false;
            mat.dg.params.full2D_solve = true;
            mat.dg.params.full2D_solver = 'bicgstab';
            mat.dg.params.full2D_preconditioner = 'blockjacobi';
            mat.dg.params.full2D_solverTol = 1e-10;
            mat.dg.params.full2D_solverMaxIt = 2000;

            DG = solve_transport_DG_full2D(mat, mat.V, 0.2, 0.1);

            testCase.verifyTrue(DG.info.solve.preconditioner.enabled);
            testCase.verifyEqual(DG.info.solve.preconditioner.requested, ...
                'blockjacobi');
            testCase.verifyEqual(DG.info.solve.preconditioner.method, ...
                'rowabs');
            testCase.verifyTrue(isfield(DG.info.solve.preconditioner, ...
                'retryFrom'));
            testCase.verifyEqual(DG.info.solve.preconditioner.retryFrom.method, ...
                'blockjacobi-rho');
            testCase.verifyEqual(DG.info.solve.flag, 0);
        end

        function autoPreconditionerRespectsStorageLimit(testCase)
            mat = makeSmallMat();
            mat.dg.params.full2D_assembleDiffMatrix = false;
            mat.dg.params.full2D_assembleDriftMatrix = false;
            mat.dg.params.full2D_solve = true;
            mat.dg.params.full2D_solver = 'bicgstab';
            mat.dg.params.full2D_preconditioner = 'auto';
            mat.dg.params.full2D_maxStoredPreconditionerDof = 100;
            mat.dg.params.full2D_solverTol = 1e-6;
            mat.dg.params.full2D_solverMaxIt = 1;

            DG = solve_transport_DG_full2D(mat, mat.V, 0.2, 0.1);

            testCase.verifyFalse(DG.info.solve.preconditioner.enabled);
            testCase.verifyEqual(DG.info.solve.preconditioner.method, 'none');
            testCase.verifyTrue(contains(DG.info.solve.preconditioner.reason, ...
                'full2D_maxStoredPreconditionerDof'));
        end

        function directSolverDoesNotAssembleAboveLimit(testCase)
            mat = makeSmallMat();
            mat.dg.params.full2D_assembleDiffMatrix = 'auto';
            mat.dg.params.full2D_assembleDriftMatrix = 'auto';
            mat.dg.params.full2D_maxAssembledDof = 100;
            mat.dg.params.full2D_solve = true;
            mat.dg.params.full2D_solver = 'direct';
            mat.dg.params.full2D_maxDirectSolveDof = 100;

            DG = solve_transport_DG_full2D(mat, mat.V, 0.2, 0.1);

            testCase.verifyEqual(DG.status, ...
                'operator-ready-direct-solve-skipped-size');
            testCase.verifyEmpty(DG.A);
            testCase.verifyEmpty(DG.rhs);
            testCase.verifyEmpty(DG.rho);
        end
    end
end

function mat = makeSmallMat()
mat = struct;
mat.type = 'Full2DDriftSystemSolveTest';
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
[X, Y] = ndgrid(mat.x, mat.y);
mat.V = 0.1*X + 0.03*Y.^2;
mat.Nd = ones(mat.Nx, mat.Ny);
mat.me_x = 0.041*ones(1, mat.Nx, mat.Ny);
mat.me_y = 0.052*ones(1, mat.Nx, mat.Ny);
mat.dg.params = struct;
mat.dg.params.full2D_coordinate_scale = 1;
mat.dg.params.N_K_chi = 3;
mat.dg.params.N_K_X = 3;
mat.dg.params.N_K_Y = 3;
mat.dg.params.N_rho_x = 5;
mat.dg.params.N_rho_y = 5;
mat.dg.params.L_rho_x = 1;
mat.dg.params.L_rho_y = 1;
mat.dg.params.full2D_fluxType = 'rusanov';
mat.dg.params.full2D_thetaLF = 0.03;
mat.dg.params.full2D_Q_diff_x = 1;
mat.dg.params.full2D_Q_diff_y = 0.7;
mat.dg.params.full2D_Q_drift = 1;
mat.dg.params.full2D_drift_scale = 1;
mat.dg.params.full2D_cap_scale = 0.1;
end

function u = deterministicVector(n)
u = reshape(sin((1:n).') + 1i*cos((1:n).'), [], 1);
end

function verifyRelativeSmall(testCase, residual, reference, tolerance)
scale = max(norm(reference), eps);
testCase.verifyLessThan(norm(residual)/scale, tolerance);
end
