function [ res ] = u_d_v_integral(grid, params,tria_index,...
    local_vertex_index, qdeg)
%U_D_V_INTEGRAL Summary of this function goes here
%   Detailed explanation goes here

f=@(llcoord) u_d_v(llocal2local(grid,local_vertex_index,llcoord), grid,...
    params,tria_index, local_vertex_index);

res = intervalquadrature(qdeg, f)*grid.EL(tria_index,local_vertex_index);

end

