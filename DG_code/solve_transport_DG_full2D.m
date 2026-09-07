function DG = solve_transport_DG_full2D(mat, Vxy, EfL, EfR)
%SOLVE_TRANSPORT_DG_FULL2D Boundary-only scaffold for full 2D Wigner-DG.
%
% Coordinates:
%   X, Y          center-of-mass coordinates, rectangular DG
%   rho_x, rho_y relative coordinates, finite volumes
%
% This function deliberately does not solve the full 4D Wigner equation yet.
% It prepares and returns the implemented boundary data/operators so they can
% be tested without assembling the large global system matrix.

addpath DG_code/core/

p = initParams_full2D(mat, Vxy, EfL, EfR);
[A, rhs, sysInfo] = get_SysM_full2D(mat, p, Vxy, EfL, EfR);

DG = struct;
DG.status = 'boundary-only';
DG.basis = 'rho-full2D';
DG.coordinates = p.index.order;
DG.centerDiscretization = p.centerDiscretization;
DG.relativeDiscretization = p.relativeDiscretization;
DG.p = p;
DG.A = A;
DG.rhs = rhs;
DG.boundary = sysInfo.boundaryData;
DG.info = sysInfo;
end
