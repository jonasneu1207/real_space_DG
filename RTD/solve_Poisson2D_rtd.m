function [Vbi, iter_err] = solve_Poisson(M, node_type, node_num, n, mat, Vbi, Vg_bias)
% SOLVE_POISSON  One step of a 2D RTD Poisson solve.
%
% Inputs:
%   M, node_type, node_num: returned by eval_Poisson(mat)
%   n      : electron density [1/m^3], size [Nx, Ny]
%   mat.N_d: donor density [1/cm^3], size [Nx, Ny]
%   mat.N_a: acceptor density [1/cm^3] (optional)
%   Vbi    : previous potential [V], size [Nx, Ny]
%   Vg_bias: right-contact boundary value [V]
%
% Contact values by default: left x contact = 0 V,
%                            right x contact = Vg_bias.
% For arbitrary Dirichlet values, supply mat.Vbc [Nx, Ny] in volts.
% The positions of Dirichlet nodes are defined by mat.boundary,
% or (if missing) automatically by eval_Poisson at the two x contacts.
%
% The sign of the source term follows the ORIGINAL code:
%     div(eps0*eps_r*grad(Vbi)) = q*(N_d - N_a - n)
% If Vbi is instead the conventional electrostatic potential phi,
% reverse the charge-source sign in BOTH branches below.
%
% Modes:
%   'direct'         : linear Poisson solve at fixed n
%   'newton-rhapson' : original spelling, Boltzmann-response Jacobian
%   'newton-raphson' : corrected spelling, same calculation
%
% iter_err = Euclidean norm of change at free nodes [V].

    q  = 1.602176634e-19;    % C
    kB = 1.380649e-23;       % J/K

    Nx = mat.Nx;
    Ny = mat.Ny;
    N  = Nx*Ny;

    assert(isequal(size(M), [N, N]), 'M must have size [Nx*Ny, Nx*Ny].');
    assert(isequal(size(n), [Nx, Ny]), 'n must have size [Nx, Ny].');
    assert(isequal(size(Vbi), [Nx, Ny]), 'Vbi must have size [Nx, Ny].');
    assert(all(isfinite(n(:))) && all(n(:) >= 0), ...
        'n must contain nonnegative, finite electron densities in m^-3.');
    assert(isfinite(Vg_bias) && isscalar(Vg_bias), ...
        'Vg_bias must be a finite scalar [V].');

    if isfield(mat, 'N_d')
        Nd = mat.N_d;
    elseif isfield(mat, 'Nd')
        Nd = mat.Nd;  % backwards compatibility
    else
        error('solve_Poisson:MissingDoping', 'mat.N_d is required.');
    end
    assert(isequal(size(Nd), [Nx,Ny]), 'mat.N_d must be [Nx, Ny].');

    if isfield(mat, 'N_a')
        Na = mat.N_a;
        assert(isequal(size(Na), [Nx,Ny]), 'mat.N_a must be [Nx, Ny].');
    else
        Na = zeros(Nx,Ny);
    end

    assert(all(isfinite(Nd(:))) && all(isfinite(Na(:))), ...
        'Doping arrays must be finite.');
    % Convert doping from cm^-3 to m^-3; n is already in m^-3.
    Nnet = (Nd - Na) * 1e6;

    % Both padded (old API) and nonpadded node arrays are accepted.
    if isequal(size(node_num), [Nx+2, Ny+2]) && ...
       isequal(size(node_type), [Nx+2, Ny+2])
        ids   = node_num(2:Nx+1, 2:Ny+1);
        types = node_type(2:Nx+1, 2:Ny+1);
    elseif isequal(size(node_num), [Nx, Ny]) && ...
           isequal(size(node_type), [Nx, Ny])
        ids   = node_num;
        types = node_type;
    else
        error('node_num and node_type have incompatible dimensions.');
    end
    assert(isequal(ids, reshape(1:N, Nx, Ny)), ...
        ['Node numbering must follow MATLAB column-major order. ', ...
         'Recompute [M,node_num,node_type] = eval_Poisson(mat).']);

    free_nodes = ids(types == 1);
    bnd_nodes  = ids(types == 2);
    assert(numel(free_nodes)+numel(bnd_nodes) == N && ...
           ~isempty(bnd_nodes), 'Invalid/free-floating node types.');

    % Dirichlet values: applied voltage on right contact, left at ground.
    if isfield(mat, 'Vbc') && ~isempty(mat.Vbc)
        Vbc = mat.Vbc;
        assert(isequal(size(Vbc), [Nx, Ny]), ...
            'mat.Vbc must have size [Nx, Ny].');
        assert(all(isfinite(Vbc(:))), 'mat.Vbc must contain finite values.');
    else
        Vbc = zeros(Nx, Ny);
        Vbc(end,:) = Vg_bias;
    end

    % MATLAB (:) is consistent with numbering returned by eval_Poisson.
    V_old = Vbi(:);
    bc = Vbc(:);
    charge = q*(Nnet(:) - n(:));

    if isfield(mat, 'poisson') && isfield(mat.poisson, 'opt')
        opt = mat.poisson.opt;
    else
        opt = struct();
    end
    if isfield(opt, 'alpha') && ~isempty(opt.alpha)
        alpha = opt.alpha;
    else
        alpha = 1;
    end
    assert(isscalar(alpha) && isfinite(alpha) && alpha > 0 && alpha <= 1, ...
        'mat.poisson.opt.alpha must lie in (0,1].');

    if isfield(opt, 'solve') && ~isempty(opt.solve)
        mode = lower(char(opt.solve));
    else
        mode = 'direct';
    end

    V_new = V_old;
    V_new(bnd_nodes) = bc(bnd_nodes);  % enforce contact potentials exactly

    switch mode
        case 'direct'
            % M_ff*V_f = charge_f - M_fb*V_b.
            rhs = charge(free_nodes) ...
                - M(free_nodes,bnd_nodes)*V_new(bnd_nodes);
            V_target = M(free_nodes,free_nodes) \ rhs;
            V_new(free_nodes) = V_old(free_nodes) ...
                + alpha*(V_target - V_old(free_nodes));

        case {'newton-rhapson', 'newton-raphson'}
            % Original Jacobian assumes dn/dVbi = -q*n/(kB*T),
            % i.e. Boltzmann electron response. This is NOT generally
            % the derivative of a quantum-transport carrier density.
            assert(isfield(mat, 'Temp') && isfinite(mat.Temp) && ...
                   isscalar(mat.Temp) && mat.Temp > 0, ...
                'mat.Temp must be positive for Newton-Raphson.');

            V0 = V_old;
            V0(bnd_nodes) = bc(bnd_nodes);
            F = M*V0 - charge;
            nf = numel(free_nodes);
            Jff = M(free_nodes,free_nodes) ...
                - spdiags(q^2*n(free_nodes)/(kB*mat.Temp), 0, nf, nf);
            dVf = -Jff \ F(free_nodes);
            V_new(free_nodes) = V0(free_nodes) + alpha*dVf;

        otherwise
            error('solve_Poisson:UnknownMethod', ...
                'Unknown Poisson solver: %s', mode);
    end

    iter_err = norm(V_new(free_nodes) - V_old(free_nodes), 2);
    Vbi = reshape(V_new, Nx, Ny);
end
