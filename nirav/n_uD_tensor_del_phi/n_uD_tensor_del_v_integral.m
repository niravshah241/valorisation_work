function [ res ] = n_uD_tensor_del_v_integral( params, grid,...
    tria_index, local_vertex_index,qdeg )

f = @(llcoord) n_uD_tensor_del_v_local (llcoord,params, grid,...
    tria_index, local_vertex_index);

res = intervalquadrature(qdeg,f) * grid.EL(tria_index,local_vertex_index);