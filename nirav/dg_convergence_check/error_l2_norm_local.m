function [ res ] = error_l2_norm_local( lcoord, params, grid, tria_index )
%ERROR_L2_NORM_LOCAL Summary of this function goes here
%   Detailed explanation goes here

basis = ldg_evaluate_basis(lcoord,params);
gids = ldg_global_dof_index(params,grid);
dofs = basis' * params.dofs(gids(tria_index,:));
glob = local2global(grid,tria_index,lcoord,params);
value = params.dof_analytical(glob);
res = norm(value-dofs,2)^2;

end