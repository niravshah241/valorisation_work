function [ reduced_basis_matrix_B, eigen_values] = ...
    pod_extension( snapshot_matrix, params, min_eigen_value, ...
    max_reduced_basis, inner_product_matrix )
%POD_EXTENSION Summary of this function goes here
%   Detailed explanation goes here

ns = size(snapshot_matrix,2);
M = inner_product_matrix;

[ V, D] = eig(snapshot_matrix' * M * snapshot_matrix);
[d,ind] = sort(diag(D),'descend');
Ds = D(ind,ind);
eigen_values = diag(Ds);
Vs = V(:,ind);

if max_reduced_basis > sum(d >= min_eigen_value)
    max_reduced_basis = sum(d >= min_eigen_value);
end

Ds = Ds(1:max_reduced_basis,1:max_reduced_basis);
Vs = Vs(:,1:max_reduced_basis);
R = [eye(max_reduced_basis);...
    zeros(ns - max_reduced_basis,max_reduced_basis)];

reduced_basis_matrix_B = zeros(params.ndofs,max_reduced_basis);
for i = 1:1:max_reduced_basis
    reduced_basis_matrix_B(:,i) = snapshot_matrix * V(:,i) ...
        / sqrt((snapshot_matrix * V(:,i))' * M * ...
        (snapshot_matrix * V(:,i)));
end

end