function results = solve(mat, solver, simulation)

switch simulation
    case 'flatband'
        results = solve2D_rtd_stat_flatband(solver, mat);

    case 'self-consistent'
        error('RTD:PeriodicPoissonNotImplemented', ...
            'The periodic-Y RTD path currently supports prescribed potentials only. Use flatband or run_RTD_full2D.');

    otherwise
        error('RTD:UnsupportedSimulation','Unknown simulation type: %s.',simulation);

end

end

