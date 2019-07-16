function [ res ] = q_n_u_d_integral( params, paramsP,grid,tria_index,...
    local_vertex_index, qdeg)
%Q_N_U_D_INTERGRAL Summary of this function goes here
%   Detailed explanation goes here

f=@(llcoord) q_n_u_d(llocal2local(grid,local_vertex_index,llcoord),params, paramsP,grid, tria_index, local_vertex_index);

res=intervalquadrature(qdeg,f)*grid.EL(tria_index,local_vertex_index);

end