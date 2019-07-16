function [ res ] = del_u_average_n_v_tensor_jump_minus_minus_integral_22...
    (params,grid,tria_index,local_vertex_index,qdeg)

if nargin == 4
    qdeg = params.qdeg;
end

f = @(llcoord) del_u_average_n_v_tensor_jump_minus_minus_local_22( llcoord,...
    params,grid,tria_index,local_vertex_index);

res = intervalquadrature(qdeg,f) * grid.EL(tria_index,local_vertex_index);

end