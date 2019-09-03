function [ res ] = error_l2_norm_integral( params, grid, tria_index, qdeg )
%ERROR_L2_NORM_INTEGRAL Summary of this function goes here
%   Detailed explanation goes here

f = @(lcoord) error_l2_norm_local(lcoord,params,grid,tria_index);
res = triaquadrature(qdeg,f)*2*grid.A(tria_index);

end