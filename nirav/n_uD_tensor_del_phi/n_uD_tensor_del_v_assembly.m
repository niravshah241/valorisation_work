function [ res ] = n_uD_tensor_del_v_assembly( params, grid, qdeg)

res = zeros(params.ndofs,1);

[ tria_index,local_edge_index,a_index, b_index, local_vertex_index] = ...
    tria_edge_index_dirichlet( grid );

gids = ldg_global_dof_index(params, grid);

for i = 1:1:length(tria_index)
    ids = gids(tria_index(i),:);
    res(ids) = res(ids) + n_uD_tensor_del_v_integral( params, grid,...
    tria_index(i), local_vertex_index(i),qdeg );
end