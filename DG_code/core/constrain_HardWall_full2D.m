function [A, rhs, info] = constrain_HardWall_full2D(A, rhs, base, p)
%CONSTRAIN_HARDWALL_FULL2D Impose homogeneous two-point Dirichlet constraints.
%
% For active injection E and exact DG mass M, the active equations are
%   (E'*M*E)^(-1)*E'*M*A*E*v = (E'*M*E)^(-1)*E'*M*b.
% The public vector size/order is retained; inactive equations are gamma*u=0.
% Those rows are algebraic constraints, not a physical absorption term.
% Apply this once AFTER adding every contribution, including potential/CAP.

wall = get_HardWallGeometry_full2D(p);
if isfield(base,'transportScales')
    scales = base.transportScales;
else
    scales = base.diff.transportScales;
end
wall.constraintScale = max(1, ...
    abs(scales.Qx)/(min(diff(p.dg.X.vertices))*p.relative.rhoX.delta) ...
    + abs(scales.Qy)/(min(diff(p.dg.Y.vertices))*p.relative.rhoY.delta));

if ~isempty(A)
    A = constrainedMatrix(A,wall);
end
if ~isempty(rhs)
    rhs = applyTransverse(wall.testProjection,rhs,wall.arraySize);
end
info = base;
info.hardWall = wall;
info.hardWallClosureDeferred = false;
info.apply = @(u) constrainedApply(base.apply,wall,u);
info.makeGPUApply = @() makeGPUApply(base,wall);
info.assemble = @() constrainedMatrix(base.assemble(),wall);
info.assembleRhs = @() applyTransverse(wall.testProjection, ...
    base.assembleRhs(),wall.arraySize);
info.getRowAbsSum = @() rowAbsBound(base,wall);
info.getDiagonal = @() constrainedDiagonal(base,wall,p);
info.getRelativeBlockData = @() constrainedBlocks(base,wall,p);
end

function y = constrainedApply(applyBase,wall,u)
activeInput = maskInput(u,wall);
y = applyTransverse(wall.testProjection,applyBase(activeInput),wall.arraySize);
y = y+wall.constraintScale*(u(:)-activeInput);
end

function v = maskInput(u,wall)
v = reshape(u,wall.arraySize).*reshape(wall.active, ...
    [1,wall.arraySize(2),1,wall.arraySize(4)]);
v = v(:);
end

function y = applyTransverse(T,u,sz)
U = permute(reshape(u,sz),[2,4,1,3]);
Y = T*reshape(U,sz(2)*sz(4),[]);
Y = ipermute(reshape(Y,sz([2,4,1,3])),[2,4,1,3]);
y = Y(:);
end

function applyGPU = makeGPUApply(base,wall)
applyBase = base.makeGPUApply();
wall.active = gpuArray(wall.active);
wall.testProjection = gpuArray(wall.testProjection);
applyGPU = @(u) constrainedApply(applyBase,wall,u);
end

function A = constrainedMatrix(A,wall)
n = prod(wall.arraySize);
active = maskInput(ones(n,1),wall);
P = spdiags(active,0,n,n);
T = liftTransverse(wall.testProjection,wall.arraySize);
A = sparse(T*A*P+wall.constraintScale*spdiags(1-active,0,n,n));
end

function T = liftTransverse(smallT,sz)
n = prod(sz);
order = permute(reshape(1:n,sz),[2,4,1,3]);
order = order(:);
[rows,cols,values] = find(kron(speye(sz(1)*sz(3)),smallT));
T = sparse(order(rows),order(cols),values,n,n);
end

function bound = rowAbsBound(base,wall)
bound = applyTransverse(abs(wall.testProjection), ...
    base.getRowAbsSum(),wall.arraySize);
active = maskInput(ones(size(bound)),wall);
bound = bound+wall.constraintScale*(1-active);
end

function blocks = constrainedBlocks(base,wall,p)
raw = base.getRelativeBlockData();
blocks = raw;
blocks.blockType = 'hard-wall-constrained-relative-blocks';
blocks.getBlock = @(id) projectedBlock(raw,wall,p,id,id);
blocks.getCoupledBlock = @(rowId,colId) projectedBlock(raw,wall,p,rowId,colId);
end

function diagonal = constrainedDiagonal(base,wall,p)
blocks = constrainedBlocks(base,wall,p);
diagonal = complex(zeros(p.relative.nDof,p.dg.nCenterDof));
for id = 1:p.dg.nCenterDof
    diagonal(:,id) = full(diag(blocks.getBlock(id)));
end
diagonal = diagonal(:);
end

function block = projectedBlock(raw,wall,p,rowId,colId)
[ix,iy] = ind2sub(p.dg.centerSize,rowId);
[~,columnY] = ind2sub(p.dg.centerSize,colId);
nR = p.relative.NrhoY;
rowIndices = (1:nR)+nR*(iy-1);
yIds = (ceil(iy/p.dg.Y.nLocal)-1)*p.dg.Y.nLocal+(1:p.dg.Y.nLocal);
block = sparse(p.relative.nDof,p.relative.nDof);
for jy = yIds
    colIndices = (1:nR)+nR*(jy-1);
    weights = full(diag(wall.testProjection(rowIndices,colIndices)));
    if any(weights)
        centerId = sub2ind(p.dg.centerSize,ix,jy);
        D = spdiags(repelem(weights,p.relative.NrhoX),0, ...
            p.relative.nDof,p.relative.nDof);
        block = block+D*raw.getCoupledBlock(centerId,colId);
    end
end
keep = repelem(double(wall.active(:,columnY)),p.relative.NrhoX);
block = block*spdiags(keep,0,p.relative.nDof,p.relative.nDof);
if rowId == colId
    block = block+wall.constraintScale*spdiags(1-keep,0, ...
        p.relative.nDof,p.relative.nDof);
end
end
