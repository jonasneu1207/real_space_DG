function [density,info] = project_DG_density_RTD(nDG,p,mat,periodicY)
%PROJECT_DG_DENSITY_RTD Conservative DG polynomial -> Poisson dual volumes.
% Uses raw nDG on BOTH DG axes, not the plotting/interpolation density n.
% Each discontinuous element is integrated separately. Periodic endpoint
% half volumes are merged; the output includes a duplicate last Y column.
arguments
    nDG {mustBeFloat,mustBeReal,mustBeFinite}
    p (1,1) struct
    mat (1,1) struct
    periodicY (1,1) logical
end
validateattributes(nDG,{'double','single'},{'size',[p.dg.X.nDof,p.dg.Y.nDof]});
if abs(p.coordinateScale-1e-9)>1e-20
    error('RTD:PoissonCoordinateUnits','The RTD Poisson material coordinates must be in nm.');
end
[Tx,wx] = integrationMap(p.dg.X,mat.x(:)*1e-9);
[Ty,wy] = integrationMap(p.dg.Y,mat.y(:)*1e-9);
if periodicY
    Ty(1,:) = Ty(1,:)+Ty(end,:);
    wy(1) = wy(1)+wy(end);
    Ty(end,:) = []; wy(end) = [];
end
volume = wx*wy.';
integral = Tx*double(nDG)*Ty.';
uniqueDensity = integral./volume;
sourceX = reshape(sum(p.dg.X.M,1).'*p.dg.X.jacobian,[],1);
sourceY = reshape(sum(p.dg.Y.M,1).'*p.dg.Y.jacobian,[],1);
sourceIntegral = sourceX.'*double(nDG)*sourceY;
targetIntegral = sum(integral(:));
density = uniqueDensity;
if periodicY, density = [density,density(:,1)]; end
info = struct('method','exact elementwise DG polynomial integrals over Poisson dual volumes', ...
    'sourceIntegral',sourceIntegral,'targetIntegral',targetIntegral, ...
    'relativeConservationError',abs(targetIntegral-sourceIntegral)/max(abs(sourceIntegral),realmin), ...
    'uniqueVolume',volume,'uniqueYNodes',numel(wy),'minimumDensity',min(density(:)), ...
    'units','density m^-3; integral per unit length in z, m^-1');
end

function [T,widths] = integrationMap(axis,nodes)
if numel(nodes)<2 || any(diff(nodes)<=0) ...
        || max(abs([nodes(1)-axis.vertices(1),nodes(end)-axis.vertices(end)]))>1e-18
    error('RTD:ChargeGridMismatch','DG and Poisson grids must cover the same interval.');
end
edges = [nodes(1);(nodes(1:end-1)+nodes(2:end))/2;nodes(end)];
widths = diff(edges);
T = sparse(numel(nodes),axis.nDof);
% n-point Gauss-Legendre integrates degree 2*n-1, including each Lagrange basis.
nq = ceil(axis.nLocal/2);
off = (1:nq-1)./sqrt(4*(1:nq-1).^2-1);
[vectors,values] = eig(diag(off,1)+diag(off,-1));
[points,order] = sort(diag(values));
weights = 2*vectors(1,order).^2;
for element = 1:axis.nElements
    left = axis.vertices(element); right = axis.vertices(element+1);
    rows = find(edges(1:end-1)<right & edges(2:end)>left);
    cols = (element-1)*axis.nLocal+(1:axis.nLocal);
    for row = rows(:).'
        a = max(edges(row),left); b = min(edges(row+1),right);
        sample = ((a+b)/2+(b-a)/2*points-(left+right)/2)/axis.jacobian(element);
        basis = ones(nq,axis.nLocal);
        for j = 1:axis.nLocal
            for k = [1:j-1,j+1:axis.nLocal]
                basis(:,j) = basis(:,j).*(sample-axis.r(k))/(axis.r(j)-axis.r(k));
            end
        end
        T(row,cols) = (b-a)/2*(weights*basis); %#ok<SPRIX> bounded local stencil, built only for the small charge grid
    end
end
end
