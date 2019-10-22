function [ res ] = ldg_h1_seminorm_mass_matrix_integral...
    ( params, grid, tria_index, qdeg )
%LDG_H1_SEMINORM_MASS_MATRIX_INTEGRAL Summary of this function goes here
%   Detailed explanation goes here

if nargin == 3
    qdeg = params.qdeg;
end

f = @(lcoord) ldg_h1_seminorm_mass_matrix_local...
    (lcoord,params,grid,tria_index);
res = triaquadrature(qdeg,f) * 2 * grid.A(tria_index);

end