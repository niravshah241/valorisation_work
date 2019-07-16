function [res] = source_local_extension(lcoord,params,paramsP,grid,k)

velocity_basis = ldg_evaluate_basis( lcoord, params);
glob = local2global( grid, k, lcoord, params);
f = params.rhs_func( glob, params, paramsP, grid);
res = velocity_basis * f;

end