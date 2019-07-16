function [ res ] = t_v( lcoord, params, paramsP, grid, ...
    tria_index, local_vertex_index)
%T_V Summary of this function goes here
%   Detailed explanation goes here

t = neumann_values(local2global(grid,tria_index,lcoord,params), ...
    params, paramsP);
v = ldg_evaluate_basis(lcoord,params);
res = v * t; 

end

