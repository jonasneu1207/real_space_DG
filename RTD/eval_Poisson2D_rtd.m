function [M, node_num, node_type] = eval_Poisson(mat)
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
% node_type: 0 = outside, 1 = free, 2 = Dirichlet.

    eps0 = 8.8541878128e-12; % F/m
    Nx = mat.Nx;
    Ny = mat.Ny;
    dx = mat.dx;              % nm
    dy = mat.dy;              % nm

    if isfield(mat, 'eps')
        er = mat.eps;
    elseif isfield(mat, 'Eps')
        er = mat.Eps;          % backwards compatibility
    else
        error('eval_Poisson:MissingEps', ...
              'mat.eps (relative permittivity) is required.');
    end

    assert(isequal(size(er), [Nx, Ny]), ...
        'mat.eps must have size [mat.Nx, mat.Ny].');
    assert(all(isfinite(er(:))) && all(er(:) > 0), ...
        'Relative permittivity must be positive and finite.');
    assert(isscalar(dx) && dx > 0 && isfinite(dx) && ...
           isscalar(dy) && dy > 0 && isfinite(dy), ...
        'mat.dx and mat.dy must be positive finite scalars (nm).');
    assert(Nx >= 2 && Ny >= 1 && Nx == floor(Nx) && Ny == floor(Ny), ...
        'Require integer Nx >= 2 and Ny >= 1.');

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
        boundary([1, Nx], :) = true;
    end

    assert(any(boundary(:)), ...
        'At least one Dirichlet node is needed to fix the potential.');

    N = Nx * Ny;
    ids = reshape(1:N, Nx, Ny);

    node_num  = zeros(Nx+2, Ny+2);
    node_type = zeros(Nx+2, Ny+2);
    node_num(2:Nx+1, 2:Ny+1) = ids;
    node_type(2:Nx+1, 2:Ny+1) = 1 + double(boundary);

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
                nz = nz + 1;
                row(nz) = k;
                col(nz) = ids(ix+1,iy);
                val(nz) = c;
                diagonal = diagonal - c;
            end

            % y-minus neighbour; if absent: zero normal flux (Neumann)
            if iy > 1
                e1 = er(ix,iy-1);
                c = cy * (2*e0*e1/(e0+e1));
                nz = nz + 1;
                row(nz) = k;
                col(nz) = ids(ix,iy-1);
                val(nz) = c;
                diagonal = diagonal - c;
            end

            % y-plus neighbour; if absent: zero normal flux (Neumann)
            if iy < Ny
                e1 = er(ix,iy+1);
                c = cy * (2*e0*e1/(e0+e1));
                nz = nz + 1;
                row(nz) = k;
                col(nz) = ids(ix,iy+1);
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
end
