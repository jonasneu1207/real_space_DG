function DG = solve_transport_DG_full2D_transient_flatband(mat, INIT)
%SOLVE_TRANSPORT_DG_FULL2D_TRANSIENT_FLATBAND Transient flatband scaffold.

DG = struct;
DG.basis = 'rho-full2D';
DG.transientMode = 'flatband';

error('DG:Full2D:NotImplemented', ...
    ['Transient flatband full 2D DG/FV transport (%s, %s) is not ', ...
     'implemented yet for grid [%d, %d] and INIT class %s.'], ...
     DG.basis, DG.transientMode, mat.Nx, mat.Ny, class(INIT));
end
