function [M, node_num, node_type, info] = eval_Poisson2D_rtd(mat)
% EVAL_POISSON  Assemble the 2D Poisson operator for an RTD.
%
% mat.eps : relative permittivity [Nx, Ny]
% mat.dx, mat.dy : grid spacing [nm]
% mat.boundary (optional): [Nx, Ny] mask, 1 = Dirichlet, 0 = free
%   Default: Dirichlet at x = first/last grid point (contacts),
%            homogeneous Neumann at y = first/last grid point.
%
% M uses MATLAB's column-major node ordering: k = ix + (iy-1)*Nx.
% Its free-node rows approximate div(eps0*eps_r*grad(V)) in SI units.
% At dielectric interfaces, harmonic means enforce flux continuity.
% Dirichlet-node rows contain the identity.
%
% node_num/node_type are (Nx+2)-by-(Ny+2) arrays, with an empty
% ghost-cell frame for compatibility with the original interface.
% mat.poisson.yBoundary='periodic': identify the two INCLUDED Y endpoints.
% Default follows full2D_Y_boundary='periodic', otherwise 'neumann'.
% node_type: 0 = outside, 1 = free, 2 = Dirichlet, 3 = periodic duplicate.
% Duplicate rows impose phi(:,end)=phi(:,1), not a second charge equation.
% mat.poisson.xBoundary='neumann': zero external X flux, including the
% endpoint half control volumes. No potential pinning; use the coupled
% Newton step with FIXED reservoir Fermi levels. Default is 'dirichlet'.
    arguments
        mat (1,1) struct
    end

    eps0 = 8.8541878128e-12; % F/m
    Nx = mat.Nx;
    Ny = mat.Ny;
    dx = mat.dx;              % nm
    dy = mat.dy;              % nm
    validateattributes(Nx,{'numeric'},{'real','finite','scalar','integer','>=',2});
    validateattributes(Ny,{'numeric'},{'real','finite','scalar','integer','positive'});

    if isfield(mat, 'eps')
        er = mat.eps;
    elseif isfield(mat, 'Eps')
        er = mat.Eps;          % backwards compatibility
    else
        error('eval_Poisson:MissingEps', ...
              'mat.eps (relative permittivity) is required.');
    end

    validateattributes(er,{'double','single'},{'real','finite','positive','size',[Nx,Ny]});
    validateattributes(dx,{'numeric'},{'real','finite','scalar','positive'});
    validateattributes(dy,{'numeric'},{'real','finite','scalar','positive'});
    assert(isscalar(dx) && dx > 0 && isfinite(dx) && ...
           isscalar(dy) && dy > 0 && isfinite(dy), ...
        'mat.dx and mat.dy must be positive finite scalars (nm).');
    assert(Nx >= 2 && Ny >= 1 && Nx == floor(Nx) && Ny == floor(Ny), ...
        'Require integer Nx >= 2 and Ny >= 1.');
    checkUniformGrid(mat,'x',Nx,dx);
    checkUniformGrid(mat,'y',Ny,dy);
    yBoundary = 'neumann';
    if isfield(mat,'dg') && isfield(mat.dg,'params') ...
            && isfield(mat.dg.params,'full2D_Y_boundary') ...
            && strcmpi(mat.dg.params.full2D_Y_boundary,'periodic')
        yBoundary = 'periodic';
    end
    if isfield(mat,'poisson') && isfield(mat.poisson,'yBoundary')
        yBoundary = lower(char(mat.poisson.yBoundary));
    end
    if ~any(strcmp(yBoundary,{'periodic','neumann'}))
        error('RTD:PoissonYBoundary','Use periodic or neumann for Poisson Y boundaries.');
    end
    periodic = strcmp(yBoundary,'periodic');
    if periodic && (Ny < 3 || max(abs(er(:,end)-er(:,1))) > 1e-12*max(er(:)))
        error('RTD:PoissonPeriodicMaterial','Periodic Y requires >=3 nodes and matching endpoint permittivities.');
    end

    xBoundary = get_RTD_poisson_x_boundary(mat);
    % Unless explicitly provided, contact nodes occupy x = 0 and x = Lx.
    if isfield(mat, 'boundary') && ~isempty(mat.boundary)
        boundary = mat.boundary;
        assert(isequal(size(boundary), [Nx, Ny]), ...
            'mat.boundary must have size [Nx, Ny].');
        assert(all(ismember(boundary(:), [0 1])), ...
            'mat.boundary must contain only 0 (free) and 1 (Dirichlet).');
        boundary = logical(boundary);
    else
        boundary = false(Nx, Ny);
        if strcmp(xBoundary,'dirichlet'), boundary([1, Nx], :) = true; end
    end

    if strcmp(xBoundary,'neumann')
        if any(boundary(:))
            error('RTD:PoissonBoundaryConflict','Neumann X requires an empty or all-free mat.boundary.');
        end
    else
        assert(any(boundary(:)), 'Dirichlet mode requires at least one fixed potential node.');
    end
    if periodic && ~isequal(boundary(:,1),boundary(:,end))
        error('RTD:PoissonPeriodicBoundary','Dirichlet masks must agree at the two Y endpoints.');
    end

    N = Nx * Ny;
    ids = reshape(1:N, Nx, Ny);

    node_num  = zeros(Nx+2, Ny+2);
    node_type = zeros(Nx+2, Ny+2);
    node_num(2:Nx+1, 2:Ny+1) = ids;
    node_type(2:Nx+1, 2:Ny+1) = 1 + double(boundary);
    if periodic
        node_type(2:Nx+1,Ny+1) = 3;
    end

    % Up to 4 neighbours + 1 diagonal per node.
    row = zeros(5*N, 1);
    col = zeros(5*N, 1);
    val = zeros(5*N, 1);
    nz = 0;

    % dx,dy are in nm; 1 nm^2 = 1e-18 m^2.
    cx = eps0 * 1e18 / dx^2;
    cy = eps0 * 1e18 / dy^2;

    for iy = 1:Ny
        for ix = 1:Nx
            k = ids(ix,iy);

            if periodic && iy == Ny
                row(nz+(1:2)) = k;
                col(nz+(1:2)) = [k;ids(ix,1)];
                val(nz+(1:2)) = [1;-1];
                nz = nz+2;
                continue;
            end

            if boundary(ix,iy)
                nz = nz + 1;
                row(nz) = k;
                col(nz) = k;
                val(nz) = 1;
                continue;
            end

            diagonal = 0;
            e0 = er(ix,iy);

            % x-minus neighbour (contact side)
            if ix > 1
                e1 = er(ix-1,iy);
                c = cx * (2*e0*e1/(e0+e1));
                if ix == Nx, c = 2*c; end % half control volume at a free X endpoint
                nz = nz + 1;
                row(nz) = k;
                col(nz) = ids(ix-1,iy);
                val(nz) = c;
                diagonal = diagonal - c;
            end

            % x-plus neighbour (opposite contact side)
            if ix < Nx
                e1 = er(ix+1,iy);
                c = cx * (2*e0*e1/(e0+e1));
                if ix == 1, c = 2*c; end
                nz = nz + 1;
                row(nz) = k;
                col(nz) = ids(ix+1,iy);
                val(nz) = c;
                diagonal = diagonal - c;
            end

            % y-minus neighbour; if absent: zero normal flux (Neumann)
            if iy > 1 || periodic
                jm = iy-1;
                if jm == 0, jm = Ny-1; end
                e1 = er(ix,jm);
                c = cy * (2*e0*e1/(e0+e1));
                if ~periodic && iy == Ny, c = 2*c; end
                nz = nz + 1;
                row(nz) = k;
                col(nz) = ids(ix,jm);
                val(nz) = c;
                diagonal = diagonal - c;
            end

            % y-plus neighbour; if absent: zero normal flux (Neumann)
            if iy < Ny
                jp = iy+1;
                if periodic && jp == Ny, jp = 1; end
                e1 = er(ix,jp);
                c = cy * (2*e0*e1/(e0+e1));
                if ~periodic && iy == 1, c = 2*c; end
                nz = nz + 1;
                row(nz) = k;
                col(nz) = ids(ix,jp);
                val(nz) = c;
                diagonal = diagonal - c;
            end

            nz = nz + 1;
            row(nz) = k;
            col(nz) = k;
            val(nz) = diagonal;
        end
    end

    M = sparse(row(1:nz), col(1:nz), val(1:nz), N, N);
    info = struct('xBoundary',xBoundary,'yBoundary',yBoundary,'uniqueYNodes',Ny-double(periodic), ...
        'potentialUnits','V','operator','div(epsilon*grad(phi))', ...
        'source','q*(n-Nd+Na), all densities in m^-3', ...
        'interfaceMean','harmonic','periodicEndpointDuplicated',periodic);
end

function checkUniformGrid(mat,name,n,spacing)
if isfield(mat,name)
    grid = mat.(name);
    if ~isreal(grid) || numel(grid) ~= n || any(~isfinite(grid(:))) ...
            || any(abs(diff(grid(:))-spacing) > 1e-10*spacing)
        error('RTD:PoissonUniformGrid','mat.%s must agree with its uniform spacing in nm.',name);
    end
end
end

