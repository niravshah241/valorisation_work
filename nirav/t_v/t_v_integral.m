function [ res ] = t_v_integral( params, paramsP, grid, ...
    tria_index, local_vertex_index, qdeg)
%Explanation: TODO

%face_index = local_vertex_index;

f = @(llcoord) t_v(llocal2local(grid,local_vertex_index,llcoord), ...
    params, paramsP, grid, tria_index, local_vertex_index);

res = intervalquadrature(qdeg,f) * grid.EL(tria_index,local_vertex_index);

end