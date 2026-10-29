function value = get_RTD_poisson_x_boundary(mat)
%GET_RTD_POISSON_X_BOUNDARY Shared contact boundary selection (default fixed).
value = 'dirichlet';
if isfield(mat,'poisson') && isfield(mat.poisson,'xBoundary')
    value = validatestring(mat.poisson.xBoundary,{'dirichlet','neumann'});
end
end
