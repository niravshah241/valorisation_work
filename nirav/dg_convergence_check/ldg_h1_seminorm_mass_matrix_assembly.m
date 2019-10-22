function [ res ] = ldg_h1_seminorm_mass_matrix_assembly( params, grid, qdeg)
%LDG_H1_SEMINORM_MASS_MATRIX_ASSEMBLY Summary of this function goes here
%   Detailed explanation goes here

if nargin == 2
    qdeg = params.qdeg;
end

res = sparse(zeros(params.ndofs));
gids = ldg_global_dof_index(params, grid);

for tria_index = 1:1:grid.nelements
    res(gids(tria_index,:),gids(tria_index,:)) = ...
        res(gids(tria_index,:),gids(tria_index,:)) + ...
        ldg_h1_seminorm_mass_matrix_integral( params, grid, ...
        tria_index, qdeg);
end

end