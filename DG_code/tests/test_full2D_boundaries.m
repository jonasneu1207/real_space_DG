function tests = test_full2D_boundaries
%TEST_FULL2D_BOUNDARIES Cheap smoke tests for the full 2D boundary scaffold.

tests = functiontests(localfunctions);
end

function setupOnce(~)
thisFile = mfilename('fullpath');
testDir = fileparts(thisFile);
dgDir = fileparts(testDir);
addpath(dgDir);
addpath(fullfile(dgDir, 'core'));
end

function testTensorProductDimensions(testCase)
mat = makeSmallMat(9, 7, 9, 7);
p = initParams_full2D(mat, mat.V, 0.1, 0.0);

verifyEqual(testCase, p.index.order, {'X-DG-DOF', 'Y-DG-DOF', 'rho_x', 'rho_y'});
verifyEqual(testCase, p.dg.X.nElements, 4);
verifyEqual(testCase, p.dg.Y.nElements, 3);
verifySize(testCase, p.dg.MXY, [p.dg.nLocalDof, p.dg.nLocalDof]);
verifySize(testCase, p.dg.DX, [p.dg.nLocalDof, p.dg.nLocalDof]);
verifySize(testCase, p.dg.DY, [p.dg.nLocalDof, p.dg.nLocalDof]);
verifySize(testCase, p.relative.Ax, [p.relative.nDof, p.relative.nDof]);
verifySize(testCase, p.relative.Ay, [p.relative.nDof, p.relative.nDof]);
end

function testRectangularDGFaces(testCase)
mat = makeSmallMat(5, 5, 5, 5);
p = initParams_full2D(mat, mat.V, 0.1, 0.0);

verifySize(testCase, p.dg.faces.local.left.trace, ...
    [p.dg.Y.nLocal, p.dg.nLocalDof]);
verifySize(testCase, p.dg.faces.local.right.trace, ...
    [p.dg.Y.nLocal, p.dg.nLocalDof]);
verifySize(testCase, p.dg.faces.local.bottom.trace, ...
    [p.dg.X.nLocal, p.dg.nLocalDof]);
verifySize(testCase, p.dg.faces.local.top.trace, ...
    [p.dg.X.nLocal, p.dg.nLocalDof]);
verifyEqual(testCase, numel(p.dg.faces.global.left.centerDofs), p.dg.Y.nDof);
verifyEqual(testCase, numel(p.dg.faces.global.right.centerDofs), p.dg.Y.nDof);
verifyEqual(testCase, numel(p.dg.faces.global.bottom.centerDofs), p.dg.X.nDof);
verifyEqual(testCase, numel(p.dg.faces.global.top.centerDofs), p.dg.X.nDof);
verifyFalse(testCase, isfield(p.dg.faces.global, 'corners'));
end

function testSourceInflowOnlyUsesSourceData(testCase)
[mat, p] = makeToyInflowProblem;
mat.dg.params.full2D_sourceRho = [10; 20; 30; 40];
mat.dg.params.full2D_drainRho = [100; 200; 300; 400];

[inflow, ~] = get_InflowBoundary_full2D(mat, p, 0.1, 0.0);
rhoInside = ones(4, p.dg.Y.nDof);
flux = apply_InflowFlux_full2D(inflow.source, rhoInside);

expected = repmat([2; 1; -30; -80], 1, p.dg.Y.nDof);
verifyEqual(testCase, flux, expected, 'AbsTol', 1e-12);
verifyEqual(testCase, inflow.source.normal, [-1, 0]);
verifyEqual(testCase, inflow.source.inflowComponents, 2);
end

function testDrainInflowOnlyUsesDrainData(testCase)
[mat, p] = makeToyInflowProblem;
mat.dg.params.full2D_sourceRho = [10; 20; 30; 40];
mat.dg.params.full2D_drainRho = [100; 200; 300; 400];

[inflow, ~] = get_InflowBoundary_full2D(mat, p, 0.1, 0.0);
rhoInside = ones(4, p.dg.Y.nDof);
flux = apply_InflowFlux_full2D(inflow.drain, rhoInside);

expected = repmat([-200; -200; 1; 2], 1, p.dg.Y.nDof);
verifyEqual(testCase, flux, expected, 'AbsTol', 1e-12);
verifyEqual(testCase, inflow.drain.normal, [1, 0]);
verifyEqual(testCase, inflow.drain.inflowComponents, 2);
end

function testMaterialDefaultReservoirIsNonzero(testCase)
mat = makeSmallMat(5, 5, 7, 5);
mat.Temp = 300;
mat.deg_factor = 1;
mat.me_x = 0.041*ones(1, mat.Nx, mat.Ny);
mat.me_y = 0.041*ones(1, mat.Nx, mat.Ny);
mat.V = zeros(mat.Nx, mat.Ny);
mat.V(end, :) = 0.05;

p = initParams_full2D(mat, mat.V, 0.2, 0.1);
[inflow, info] = get_InflowBoundary_full2D(mat, p, 0.2, 0.1, mat.V);

verifyEqual(testCase, info.source.origin, 'material-default');
verifyEqual(testCase, info.drain.origin, 'material-default');
verifySize(testCase, inflow.source.rhoBoundary, [p.relative.nDof, p.dg.Y.nDof]);
verifySize(testCase, inflow.drain.rhoBoundary, [p.relative.nDof, p.dg.Y.nDof]);
verifyGreaterThan(testCase, norm(inflow.source.rhoBoundary, 'fro'), 0);
verifyGreaterThan(testCase, norm(inflow.drain.rhoBoundary, 'fro'), 0);
end

function testSpecularReflectionIsInvolutionAndSparse(testCase)
mat = makeSmallMat(5, 5, 5, 7);
p = initParams_full2D(mat, mat.V, 0.1, 0.0);
reflection = get_SpecularReflection_full2D(p);

I = speye(p.relative.nDof);
verifyTrue(testCase, issparse(reflection.R_y));
verifyEqual(testCase, nnz(reflection.R_y), p.relative.nDof);
verifyLessThan(testCase, norm(reflection.R_y*reflection.R_y - I, 'fro'), 1e-12);

rho = reshape(1:p.relative.nDof, [], 1);
verifyEqual(testCase, reflection.R_y*(reflection.R_y*rho), rho);
end

function testCAPOnlyActsNearRelativeBoundariesAndIsSparse(testCase)
mat = makeSmallMat(5, 5, 21, 19);
p = initParams_full2D(mat, mat.V, 0.1, 0.0);
cap = get_CAP_full2D(p);

[maskX, maskY] = ndgrid(cap.maskX, cap.maskY);
interiorMask = ~(maskX | maskY);

verifyTrue(testCase, issparse(cap.Crho));
verifyTrue(testCase, issparse(cap.potentialMatrix));
verifyGreaterThan(testCase, nnz(cap.Crho), 0);
verifyEqual(testCase, cap.profile2D(interiorMask), ...
    zeros(nnz(interiorMask), 1), 'AbsTol', 1e-14);
verifyEqual(testCase, cap.mask2D, maskX | maskY);
end

function testBoundaryAggregatorKeepsTypesSeparate(testCase)
mat = makeSmallMat(5, 5, 7, 5);
p = initParams_full2D(mat, mat.V, 0.1, 0.0);
[boundary, info] = get_Boundary_full2D(mat, p, 0.1, 0.0);

verifyEqual(testCase, boundary.physical.XLeft.type, 'characteristic-inflow');
verifyEqual(testCase, boundary.physical.XRight.type, 'characteristic-inflow');
verifyEqual(testCase, boundary.physical.YBottom.type, 'specular-reflection');
verifyEqual(testCase, boundary.physical.YTop.type, 'specular-reflection');
verifyTrue(testCase, isfield(boundary.relative, 'CAP'));
verifyGreaterThan(testCase, info.capNnz, 0);
verifyGreaterThan(testCase, info.reflectionNnz, 0);
end

function mat = makeSmallMat(nx, ny, nRhoX, nRhoY)
mat = struct;
mat.type = 'Full2DTest';
mat.Nx = nx;
mat.Ny = ny;
mat.x = linspace(0, 1, nx);
mat.y = linspace(-0.5, 0.5, ny);
mat.dx = mat.x(2)-mat.x(1);
mat.dy = mat.y(2)-mat.y(1);
mat.V = zeros(nx, ny);
mat.dg.params = struct;
mat.dg.params.full2D_coordinate_scale = 1;
mat.dg.params.N_K_chi = 3;
mat.dg.params.N_K_X = 3;
mat.dg.params.N_K_Y = 3;
mat.dg.params.N_rho_x = nRhoX;
mat.dg.params.N_rho_y = nRhoY;
mat.dg.params.L_rho_x = 1;
mat.dg.params.L_rho_y = 1;
end

function [mat, p] = makeToyInflowProblem
mat = makeSmallMat(3, 3, 4, 3);
p = struct;
p.domain.normals.XLeft = [-1, 0];
p.domain.normals.XRight = [1, 0];
p.domain.normalConvention = ...
    'Normals point out of the physical X-Y domain; negative eigenvalues of A_n are inflow.';
p.index.order = {'X-DG-DOF', 'Y-DG-DOF', 'rho_x', 'rho_y'};
p.relative.NrhoX = 4;
p.relative.NrhoY = 1;
p.relative.size = [4, 1];
p.relative.nDof = 4;
p.relative.Ax1D = spdiags([-2; -1; 1; 2], 0, 4, 4);
p.relative.Ay1D = sparse(1, 1);
p.relative.Ax = p.relative.Ax1D;
p.relative.Ay = sparse(4, 4);
p.relative.identity = speye(4);
p.dg.Y.nDof = 3;
p.dg.faces.global.left.centerDofs = (1:3).';
p.dg.faces.global.right.centerDofs = (1:3).';
end
