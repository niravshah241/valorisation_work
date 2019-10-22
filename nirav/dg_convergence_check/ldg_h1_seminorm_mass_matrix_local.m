function [ res ] = ldg_h1_seminorm_mass_matrix_local...
    (lcoord,params,grid,tria_index)
%LDG_H1_SEMINORM_MASS_MATRIX_LOCAL Summary of this function goes here
%   Detailed explanation goes here

res = zeros(params.ndofs_per_element);
basis_derivative = ldg_evaluate_basis_derivative(lcoord,params);
for i = 1:1:params.ndofs_per_element
    for j = 1:1:params.ndofs_per_element
        res(i,j) = res(i,j) + ...
            sum(sum(basis_derivative{i}.*basis_derivative{j}));
    end
end
end