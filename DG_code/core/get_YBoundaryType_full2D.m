function boundaryType = get_YBoundaryType_full2D(mat)
%GET_YBOUNDARYTYPE_FULL2D Resolve the physical top/bottom boundary model.
%
% Supported values of mat.dg.params.full2D_Y_boundary:
%   'specular'   Incoming characteristics are the rho_y-reflected interior
%                state. This is the backwards-compatible default.
%   'zero-inflow' Incoming characteristics are zero, while outgoing
%                characteristics remain determined by the interior DG
%                solution. This is an open/absorbing diagnostic boundary,
%                not a Dirichlet-zero condition on the complete rho state.

value = 'specular';
if isstruct(mat) && isfield(mat, 'dg') && isstruct(mat.dg) ...
        && isfield(mat.dg, 'params') && isstruct(mat.dg.params)
    params = mat.dg.params;
    if isfield(params, 'full2D_Y_boundary') ...
            && ~isempty(params.full2D_Y_boundary)
        value = params.full2D_Y_boundary;
    elseif isfield(params, 'full2D') && isstruct(params.full2D) ...
            && isfield(params.full2D, 'full2D_Y_boundary') ...
            && ~isempty(params.full2D.full2D_Y_boundary)
        value = params.full2D.full2D_Y_boundary;
    end
end

boundaryType = lower(strtrim(char(value)));
if any(strcmp(boundaryType, ...
        {'specular-reflection', 'reflection', 'reflecting'}))
    boundaryType = 'specular';
elseif any(strcmp(boundaryType, ...
        {'zero_inflow', 'zeroinflow', 'zero'}))
    boundaryType = 'zero-inflow';
end

if ~any(strcmp(boundaryType, {'specular', 'zero-inflow'}))
    error('DG:Full2D:UnknownYBoundaryType', ...
        ['Unknown full2D_Y_boundary "%s". Use "specular" or ', ...
         '"zero-inflow".'], boundaryType);
end
end
