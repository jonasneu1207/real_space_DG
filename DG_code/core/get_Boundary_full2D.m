function [boundary, info] = get_Boundary_full2D(mat, p, EfL, EfR, Vxy)
%GET_BOUNDARY_FULL2D Boundary data for the full physical 2D Wigner-DG path.
%
% Implemented boundary types for the rectangular X-Y domain:
%   X-left   Source characteristic inflow, outward normal (-1,0)
%   X-right  Drain characteristic inflow, outward normal ( 1,0)
%   Y-bottom/top: specular, zero-inflow, or two-point Dirichlet hard-wall.
%
% In hard-wall mode the two-point Dirichlet constraints also take precedence
% at X/Y corners. Other modes retain the adjacent-face contributions.

if nargin < 5
    Vxy = [];
end

[inflow, inflowInfo] = get_InflowBoundary_full2D(mat, p, EfL, EfR, Vxy);
reflection = get_SpecularReflection_full2D(p);
cap = get_CAP_full2D(p);
yBoundaryType = get_YBoundaryType_full2D(mat);

boundary = struct;
boundary.order = {'X-left', 'X-right', 'Y-bottom', 'Y-top'};
boundary.physical.XLeft = inflow.source;
boundary.physical.XRight = inflow.drain;
boundary.physical.YBottom = reflection.bottom;
boundary.physical.YTop = reflection.top;
boundary.physical.YBottom.R_y = reflection.R_y;
boundary.physical.YTop.R_y = reflection.R_y;
if strcmp(yBoundaryType, 'specular')
    physicalYType = reflection.type;
    incomingState = 'rho_in = R_y*rho_inside';
elseif strcmp(yBoundaryType, 'hard-wall')
    physicalYType = 'two-point-dirichlet-hard-wall';
    incomingState = 'rho = 0 if Y +/- rho_y/2 reaches either wall';
    boundary.physical.YBottom = rmfield(boundary.physical.YBottom,{'ghostOperator','R_y'});
    boundary.physical.YTop = rmfield(boundary.physical.YTop,{'ghostOperator','R_y'});
else
    physicalYType = 'characteristic-zero-inflow';
    incomingState = 'rho_in = 0; rho_out comes from the interior trace';
end
boundary.physical.YBottom.type = physicalYType;
boundary.physical.YTop.type = physicalYType;
boundary.physical.YBottom.incomingState = incomingState;
boundary.physical.YTop.incomingState = incomingState;
boundary.yBoundaryType = yBoundaryType;

% The CAP lives on the relative rho_x/rho_y boundary and is deliberately
% separate from the physical X/Y boundary fluxes.
boundary.relative.CAP = cap;
boundary.specularReflection = reflection;

info = struct;
info.full2D = true;
info.inflow = inflowInfo;
info.capNnz = nnz(cap.Crho);
info.yBoundaryType = yBoundaryType;
info.reflectionActive = strcmp(yBoundaryType, 'specular');
info.hardWallActive = strcmp(yBoundaryType, 'hard-wall');
info.reflectionNnz = info.reflectionActive*nnz(reflection.R_y);
info.availableReflectionNnz = nnz(reflection.R_y);
info.normalConvention = p.domain.normalConvention;
info.dofOrder = p.index.order;
end
