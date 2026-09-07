function DG = solve_transport_DG_full2D_transient_SC_new(mat, INIT)
%SOLVE_TRANSPORT_DG_FULL2D_TRANSIENT_SC_NEW Self-consistent transient scaffold.

DG = struct;
DG.basis = 'rho-full2D';
DG.transientMode = 'self-consistent';

error('DG:Full2D:NotImplemented', ...
    ['Self-consistent transient full 2D DG/FV transport (%s, %s) is ', ...
     'not implemented yet for grid [%d, %d] and INIT class %s.'], ...
     DG.basis, DG.transientMode, mat.Nx, mat.Ny, class(INIT));
end
