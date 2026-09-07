function [boundary, info] = get_Boundary_full2D(mat, p, EfL, EfR, Vxy)
%GET_BOUNDARY_FULL2D Boundary data for the full physical 2D Wigner-DG path.
%
% Implemented boundary types for the rectangular X-Y domain:
%   X-left   Source characteristic inflow, outward normal (-1,0)
%   X-right  Drain characteristic inflow, outward normal ( 1,0)
%   Y-bottom Specular reflection, outward normal (0,-1)
%   Y-top    Specular reflection, outward normal (0, 1)
%
% No special corner condition is introduced. Corner DOFs receive only the
% contributions from their adjacent tensor-product faces.

if nargin < 5
    Vxy = [];
end

[inflow, inflowInfo] = get_InflowBoundary_full2D(mat, p, EfL, EfR, Vxy);
reflection = get_SpecularReflection_full2D(p);
cap = get_CAP_full2D(p);

boundary = struct;
boundary.order = {'X-left', 'X-right', 'Y-bottom', 'Y-top'};
boundary.physical.XLeft = inflow.source;
boundary.physical.XRight = inflow.drain;
boundary.physical.YBottom = reflection.bottom;
boundary.physical.YTop = reflection.top;
boundary.physical.YBottom.type = reflection.type;
boundary.physical.YTop.type = reflection.type;
boundary.physical.YBottom.R_y = reflection.R_y;
boundary.physical.YTop.R_y = reflection.R_y;

% The CAP lives on the relative rho_x/rho_y boundary and is deliberately
% separate from the physical X/Y boundary fluxes.
boundary.relative.CAP = cap;
boundary.specularReflection = reflection;

info = struct;
info.full2D = true;
info.inflow = inflowInfo;
info.capNnz = nnz(cap.Crho);
info.reflectionNnz = nnz(reflection.R_y);
info.normalConvention = p.domain.normalConvention;
info.dofOrder = p.index.order;
end
