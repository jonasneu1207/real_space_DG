function INIT = get_init_DG_full2D_trans_flatband(mat, Vg)
%GET_INIT_DG_FULL2D_TRANS_FLATBAND Flatband init scaffold for full 2D DG.

INIT = struct;
INIT.basis = 'rho-full2D';
INIT.transientMode = 'flatband';

error('DG:Full2D:NotImplemented', ...
    ['Flatband initialization for the full 2D DG/FV solver (%s, %s) ', ...
     'is not implemented yet for grid [%d, %d] and Vg=%g.'], ...
     INIT.basis, INIT.transientMode, mat.Nx, mat.Ny, Vg);
end
