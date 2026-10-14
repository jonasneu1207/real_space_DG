function results = solve2D_rtd_self_consistent(SOLVER, mat)
if isfield(mat,'dg') && isfield(mat.dg,'params') ...
        && isfield(mat.dg.params,'full2D_Y_boundary') ...
        && strcmpi(mat.dg.params.full2D_Y_boundary,'periodic')
    error('RTD:PeriodicPoissonNotImplemented', ...
        'The prepared Poisson iteration has not been adapted to periodic Y. Use run_RTD_full2D with a prescribed potential.');
end
% solveIABC_stat_flatband solves for the differenent solvers SOLVER and
% material structures mat

%% solve for the equilibrium fermi level
[EfL, EfR] = getFermiLevel(mat);

% built in bias due to doping concentration
Vbi_dop = EfL-EfR;

%% data decompression
x = mat.x;
dx = mat.dx;
Nx = length(x);
Nc1 = ceil(mat.Lc1/dx);
Nc2 = ceil(mat.Lc2/dx);

TOLER = mat.poisson.opt.iter_err;
MAXIT = mat.poisson.opt.iter_max;

V_DD = mat.Vdd;
N_DD = length(V_DD);

% flag to accelerate the wigner comutations massively -> only drift part is
% re-computed to start the flag must be 0 it turns to 1 automatically!!!!
% after the first iteration has been performed
fastened_flag = 0;

%% loop over applied biases

Id_sweep         = zeros(Nx, N_DD);
n_electron_sweep = zeros(Nx, N_DD);
Vbi_sweep        = zeros(Nx, N_DD);



for I_DD = 1 : N_DD
    
    fprintf('[BIAS SWEEP AT :\t %d of %d]\n', I_DD, length(V_DD));
    
    % prepare the guess potential or flatband potential
    % lineare interpolate from known solutions --> increased numerical
    % stability as well as shorter runtimes
    if I_DD < 3
        
        Vbi = [zeros(1, Nc1), linspace(0, 1, Nx-(Nc1+Nc2)), ones(1, Nc2)].*(Vbi_dop-V_DD(I_DD));
        
    else
        
        Vbi = 2*Vbi_sweep(:,I_DD-1) - Vbi_sweep(:,I_DD-2);
        Vbi = reshape(Vbi, 1, Nx);
    
    end    
    figure(666), plot(mat.x, Vbi), drawnow;
    % iteration error
    iter_err = 5;   % assign anything greater TOLER to ensure start
    
    % iteration counter
    iter_count = 0;
    
    while iter_err > TOLER && iter_count < MAXIT
        
        results = solve_transport_DG_full2D(mat, mat.Vc, EfL, EfR);
           
                
            
                
         % SWITCH.SOLVER
        
        % solve the POISSON EQUATION
        p_hole = zeros(size(n_electron_sweep(:, I_DD))); % no need for hole concentration in this solver!
        [Vbi, iter_err] = solve_POISSON(n_electron_sweep(:, I_DD), p_hole, Vbi, V_DD(I_DD), mat, SOLVER);
        iter_count = iter_count + 1;
        
        fprintf('.....inner poisson loop at: \t %2d with error : %1.2e\n', iter_count, iter_err);
        
    end % WHILE SELF-CONSISTENT LOOP
    
    Vbi_sweep(:, I_DD) = Vbi;
    figure(9), plotyy(mat.x, Vbi+mat.Vc, mat.x, n_electron_sweep(:,I_DD)); drawnow;
    

% global output
results.Id = Id_sweep;
results.n_electron = n_electron_sweep;
results.dist_func_sweep = dist_func_sweep;
results.Vc = Vbi_sweep+mat.Vc';

end % function

