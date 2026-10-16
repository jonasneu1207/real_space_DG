function results = solve(mat, solver, simulation, options)

if nargin < 4, options = struct; end

switch simulation
    case 'flatband'
        results = solve2D_rtd_stat_flatband(solver, mat);

    case 'self-consistent'
        results = solve2D_rtd_self_consistent(solver, mat, options);

    otherwise
        error('RTD:UnsupportedSimulation','Unknown simulation type: %s.',simulation);

end

end

