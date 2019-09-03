function [ res ] = error_rb_l2_norm_local( lcoord, params, grid, ...
    reduced_dofs, tria_index )
%ERROR_RB_L2_NORM_LOCAL Summary of this function goes here
%   Detailed explanation goes here

basis = ldg_evaluate_basis(lcoord,params);
gids = ldg_global_dof_index(params,grid);
error = basis' * ( params.dofs(gids(tria_index,:)) - ...
    reduced_dofs(gids(tria_index,:)));
res = norm(error,2)^2;

end