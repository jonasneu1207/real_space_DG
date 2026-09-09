function p = initParams_full2D(mat, Vxy, EfL, EfR)
%INITPARAMS_FULL2D Parameters for the full physical 2D Wigner-DG scaffold.
%
% Global DOF order:
%   F(iRhoX, iRhoY, iX, iY) is vectorized as F(:). Thus rho_x FV cells are
%   the fastest index, followed by rho_y FV cells, X-DG DOFs and Y-DG DOFs.
%   A separable operator O_X is lifted as
%   kron(I_Y, kron(O_X, kron(I_rho_y, I_rho_x))).
%
% Center coordinates X and Y use rectangular DG elements with tensor-product
% basis phi_ij(X,Y) = phi_i^X(X)*phi_j^Y(Y). Relative coordinates rho_x and
% rho_y use finite-volume cell centers.

p = struct;
p.centerCoordinates = {'X', 'Y'};
p.relativeCoordinates = {'rho_x', 'rho_y'};
p.centerDiscretization = 'DG';
p.relativeDiscretization = 'FV';
p.EfL = EfL;
p.EfR = EfR;
p.VxySize = size(Vxy);

params = struct;
if isfield(mat, 'dg') && isfield(mat.dg, 'params')
    params = mat.dg.params;
end

coordScale = readFull2DParam(params, 'full2D_coordinate_scale', 1e-9);
p.coordinateScale = coordScale;
x = getCoordinateVector(mat, 'x', 'dx', 'Nx')*coordScale;
y = getCoordinateVector(mat, 'y', 'dy', 'Ny')*coordScale;

p.domain.X = [x(1), x(end)];
p.domain.Y = [y(1), y(end)];
p.domain.normals.XLeft = [-1, 0];
p.domain.normals.XRight = [1, 0];
p.domain.normals.YBottom = [0, -1];
p.domain.normals.YTop = [0, 1];
p.domain.normalConvention = ...
    'Normals point out of the physical X-Y domain; negative eigenvalues of A_n are inflow.';

NKDefault = readFull2DParam(params, 'N_K_chi', 3);
N_K_X = readFull2DParam(params, 'N_K_X', NKDefault);
N_K_Y = readFull2DParam(params, 'N_K_Y', NKDefault);
p.dg.X = makeRectDGAxis(x, N_K_X, params, 'X');
p.dg.Y = makeRectDGAxis(y, N_K_Y, params, 'Y');
p.dg.nLocalDof = p.dg.X.nLocal*p.dg.Y.nLocal;
p.dg.centerSize = [p.dg.X.nDof, p.dg.Y.nDof];
p.dg.nCenterDof = prod(p.dg.centerSize);

% Local tensor-product operators on one rectangular element. With local
% ordering (X-local, Y-local), X remains the fastest index.
p.dg.MXY = kron(p.dg.Y.M, p.dg.X.M);
p.dg.DX = kron(speye(p.dg.Y.nLocal), p.dg.X.Dr);
p.dg.DY = kron(p.dg.Y.Dr, speye(p.dg.X.nLocal));
p.dg.SX = kron(p.dg.Y.M, p.dg.X.S);
p.dg.SY = kron(p.dg.Y.S, p.dg.X.M);
p.dg.faces = makeRectDGFaces(p.dg.X, p.dg.Y);

N_rho_x = readFull2DParam(params, 'N_rho_x', readFull2DParam(params, 'N_xi', 40));
N_rho_y = readFull2DParam(params, 'N_rho_y', N_rho_x);
L_rho_default = readFull2DParam(params, 'Ly', max(abs(y(end)-y(1)), eps));
L_rho_x = readFull2DParam(params, 'L_rho_x', L_rho_default);
L_rho_y = readFull2DParam(params, 'L_rho_y', L_rho_default);

p.relative.rhoX = makeFVRelativeAxis(N_rho_x, L_rho_x, 'rho_x');
p.relative.rhoY = makeFVRelativeAxis(N_rho_y, L_rho_y, 'rho_y');
p.relative.NrhoX = p.relative.rhoX.nCells;
p.relative.NrhoY = p.relative.rhoY.nCells;
p.relative.size = [p.relative.NrhoX, p.relative.NrhoY];
p.relative.nDof = prod(p.relative.size);

% In rho representation, multiplication by k_x/k_y corresponds to
% -1i*d/d(rho_x/rho_y). The 2D relative operators are Kronecker sums and
% stay sparse in the FV basis.
Ix = speye(p.relative.NrhoX);
Iy = speye(p.relative.NrhoY);
p.relative.DrhoX = p.relative.rhoX.D;
p.relative.DrhoY = p.relative.rhoY.D;
p.relative.Ax1D = -1i*p.relative.DrhoX;
p.relative.Ay1D = -1i*p.relative.DrhoY;
p.relative.Ax = kron(Iy, p.relative.Ax1D);
p.relative.Ay = kron(p.relative.Ay1D, Ix);
p.relative.identity = speye(p.relative.nDof);

p.index.order = {'rho_x', 'rho_y', 'X-DG-DOF', 'Y-DG-DOF'};
p.index.arraySize = [p.relative.size, p.dg.centerSize];
p.index.nTotal = prod(p.index.arraySize);
p.index.lift = @(OX, OY, ORX, ORY) kron(OY, kron(OX, kron(ORY, ORX)));
end

function value = readFull2DParam(params, name, defaultValue)
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

function coord = getCoordinateVector(mat, coordName, spacingName, countName)
if isfield(mat, coordName) && ~isempty(mat.(coordName))
    coord = mat.(coordName)(:).';
else
    n = mat.(countName);
    d = mat.(spacingName);
    coord = (0:n-1)*d;
end
if numel(coord) < 2
    error('DG:Full2D:InvalidGrid', '%s needs at least two grid points.', coordName);
end
end

function axis = makeRectDGAxis(coord, nLocal, params, label)
if nLocal < 2
    error('DG:Full2D:InvalidDGOrder', 'N_K_%s must be at least 2.', label);
end

nIntervals = numel(coord)-1;
nElementField = ['N_el_', label];
requestedElements = readFull2DParam(params, nElementField, []);
if ~isempty(requestedElements)
    nElements = requestedElements;
else
    if mod(nIntervals, nLocal-1) ~= 0
        error('DG:Full2D:InvalidDGGrid', ...
            ['%s grid has %d intervals. This is not compatible with %d ', ...
             'local DG nodes. Set mat.dg.params.%s explicitly.'], ...
             label, nIntervals, nLocal, nElementField);
    end
    nElements = nIntervals/(nLocal-1);
end

[r, w] = JacobiGL(0, 0, nLocal-1);
V = Vandermonde1D(nLocal-1, r);
invV = inv(V);
M = sparse(invV.'*invV);
Dr = sparse(Dmatrix1D(nLocal-1, r, V));
S = sparse(M*Dr);

vertices = linspace(coord(1), coord(end), nElements+1);
nodes = zeros(1, nElements*nLocal);
jacobian = zeros(1, nElements);
for elem = 1:nElements
    xL = vertices(elem);
    xR = vertices(elem+1);
    jacobian(elem) = (xR-xL)/2;
    ids = (elem-1)*nLocal+1:elem*nLocal;
    nodes(ids) = (xL+xR)/2 + jacobian(elem)*r.';
end

axis.label = label;
axis.nLocal = nLocal;
axis.nElements = nElements;
axis.nDof = nElements*nLocal;
axis.r = r;
axis.w = w;
axis.V = V;
axis.M = M;
axis.Dr = Dr;
axis.S = S;
axis.vertices = vertices;
axis.nodes = nodes;
axis.jacobian = jacobian;
axis.faceLeft = sparse(1, 1, 1, 1, nLocal);
axis.faceRight = sparse(1, nLocal, 1, 1, nLocal);
axis.globalMass = kron(speye(nElements), M);
end

function axis = makeFVRelativeAxis(nCells, lengthValue, label)
if nCells < 2
    error('DG:Full2D:InvalidRelativeGrid', '%s needs at least two FV cells.', label);
end
edges = linspace(-lengthValue/2, lengthValue/2, nCells+1);
cells = 0.5*(edges(1:end-1)+edges(2:end));
delta = edges(2)-edges(1);

axis.label = label;
axis.nCells = nCells;
axis.length = lengthValue;
axis.edges = edges(:);
axis.cells = cells(:);
axis.delta = delta;
axis.D = spdiags(ones(nCells, 1)*[-1/(2*delta), 1/(2*delta)], [-1, 1], nCells, nCells);
end

function faces = makeRectDGFaces(axisX, axisY)
NX = axisX.nLocal;
NY = axisY.nLocal;
IX = speye(NX);
IY = speye(NY);

% Local face traces on a rectangle. They select edge values from one
% tensor-product element without introducing a special corner condition.
faces.local.left.trace = kron(IY, axisX.faceLeft);
faces.local.right.trace = kron(IY, axisX.faceRight);
faces.local.bottom.trace = kron(axisY.faceLeft, IX);
faces.local.top.trace = kron(axisY.faceRight, IX);
faces.local.left.normal = [-1, 0];
faces.local.right.normal = [1, 0];
faces.local.bottom.normal = [0, -1];
faces.local.top.normal = [0, 1];
faces.local.left.mass = axisY.M;
faces.local.right.mass = axisY.M;
faces.local.bottom.mass = axisX.M;
faces.local.top.mass = axisX.M;

nXDof = axisX.nDof;
nYDof = axisY.nDof;
faces.global.left.centerDofs = sub2ind([nXDof, nYDof], ones(nYDof, 1), (1:nYDof).');
faces.global.right.centerDofs = sub2ind([nXDof, nYDof], nXDof*ones(nYDof, 1), (1:nYDof).');
faces.global.bottom.centerDofs = sub2ind([nXDof, nYDof], (1:nXDof).', ones(nXDof, 1));
faces.global.top.centerDofs = sub2ind([nXDof, nYDof], (1:nXDof).', nYDof*ones(nXDof, 1));
end
