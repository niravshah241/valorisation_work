function [ res ] = u_d_v( lcoord, grid, params,tria_index, local_vertex_index)
%U_D_V Summary of this function goes here
%   Detailed explanation goes here

v = ldg_evaluate_basis(lcoord,params);
u_d = dirichlet_values(local2global(grid,tria_index,lcoord,params),params);

res=v*u_d;

end