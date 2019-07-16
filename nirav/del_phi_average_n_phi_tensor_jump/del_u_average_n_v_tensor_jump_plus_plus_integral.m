function [ res ] = del_u_average_n_v_tensor_jump_plus_plus_integral...
    (params,grid,tria_index,local_vertex_index,qdeg)

f = @(llcoord) del_u_average_n_v_tensor_jump_plus_plus_local( llcoord,...
    params,grid,tria_index,local_vertex_index);

res = intervalquadrature(qdeg,f) * grid.EL(tria_index,local_vertex_index);