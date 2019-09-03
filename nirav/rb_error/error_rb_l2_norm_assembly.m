function [ res ] = error_rb_l2_norm_assembly( params, grid, ...
    reduced_dofs, qdeg )

if nargin == 3
    qdeg = params.qdeg;
end

res = zeros(grid.nelements,1);

for tria_index = 1:1:grid.nelements
    res_integral = error_rb_l2_norm_integral( params, grid, ...
        reduced_dofs, tria_index, qdeg );
    res(tria_index) = res(tria_index) + res_integral;
end

res = (sum(res))^(1/2);

end