%% Periodic-Y RTD, prescribed potential (no Poisson self-consistency)
% Run from this folder or by its absolute path. No automatic MAT-file writes.
addpath(fileparts(mfilename('fullpath')));
options = struct;
options.Ly = 3;                 % nm, physical transverse period W
options.dy = 0.25;               % nm, material-grid sampling
options.simulation = 'self-consistent';       % 'self-consistent' or 'flatband'
options.selfConsistent.poissonMethod = 'newton-raphson';     % 'newton_raphson' or 'direct'
options.selfConsistent.poissonXBoundary = 'dirichlet'; % 'dirichlet' or 'neumann'; Y stays periodic
% options.selfConsistent.chargeNeutralityTolerance = 1e-3; % Neumann total charge criterion
options.selfConsistent.mixing = 0.3;
options.bias = 0.0:0.02:0.6;               % V; a vector enables a fixed-potential sweep
options.dg.N_el_X = 120;
options.dg.N_el_Y = 2;          % even, so the two-particle half shift is exact
options.dg.N_rho_x = 121;
options.dg.N_rho_y = 12;        % even; L_rho_y = 2*W is set by the runner
options.dg.L_rho_x = 100e-9;     % m
options.dg.full2D_solver = 'gmres';
options.dg.full2D_preconditioner = 'line-transverse';
options.dg.full2D_gmresRestart = 40; % RTD starting value, not a universal optimum
options.dg.full2D_lineRelativeShift = 0.003;
options.dg.full2D_solverTol = 5e-6;
options.dg.full2D_gpu = true;  % enable manually on the GPU machine
% options.removeBarriers = true; % homogeneous equilibrium check
% options.dg.full2D_periodicModeOrder = 0; % only k_y=0, clean 1D comparison
% options.dg.full2D_cap_scale = 0; % disable also rho_x CAP for a control run

%% Test flux options for more performance
% options.dg.full2D_X_flux = 'rusanov';
% options.dg.full2D_Y_flux = 'rusanov';
% options.dg.full2D_thetaLF_Y = 1.0;
% options.dg.full2D_thetaLF_X = 1.0;
%%

[results,mat] = run_RTD_full2D(options);
fprintf('RTD: flag=%g, relres=%.3e, periodicity defect=%.3e\n', ...
    results.info.solve.flag,results.info.solve.relres, ...
    results.info.solve.periodicYStateError);

% This utility selects the nearest FV cell to rho_y=0, not an interpolation.
% [S,ax] = results.getRhoSlice2D({'X','rho_x'}, ...
%     struct('Y',options.Ly*0.5e-9,'rho_y',0));
% figure;
% mesh(ax.rho_x*1e9,ax.X*1e9,real(S));
% xlabel('\rho_x [nm]'); ylabel('X [nm]'); zlabel('Re(\rho)');
% title('RTD: periodic Y, constant effective mass');


ordnername = char(datetime('now', 'Format', 'yyyy-MM-dd_HH-mm'));

speicherpfad = fullfile(pwd, ordnername);

mkdir(speicherpfad);

save(fullfile(speicherpfad, 'Sim_Res.mat'), '-struct', 'results', 'sweep', 'selfConsistency', 'mat', 'boundary');


