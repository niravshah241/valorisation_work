function [ res ] = source_assembly_extension( params, paramsP, grid )

res = zeros(params.ndofs,1);
gids = ldg_global_dof_index(params, grid);

for k = 1:1:grid.nelements
    A = source_integral_extension(params,paramsP,grid,k);
    res(gids(k,:)) = res(gids(k,:)) + A;
end

end