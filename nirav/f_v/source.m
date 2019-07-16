function res=source(lcoord,params, paramsP, grid, k)

glob=local2global(grid,k,lcoord,params);
f=params.rhs_func(glob,params,paramsP,grid);
res=ldg_evaluate_basis(lcoord,params)*f;

end