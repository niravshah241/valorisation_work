function [ res ] = error_l2_norm_assembly( params, grid, qdeg )
%ERROR_L2_NORM_ASSEMBLY Summary of this function goes here
%   Detailed explanation goes here

if nargin == 2
    qdeg = params.qdeg;
end

res = zeros(grid.nelements,1);

for tria_index = 1:1:grid.nelements
    res_integral = error_l2_norm_integral( params, grid, tria_index, qdeg );
    res(tria_index) = res(tria_index) + res_integral;
end

res = (sum(res))^(1/2);