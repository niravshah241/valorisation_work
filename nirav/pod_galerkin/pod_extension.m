function [ reduced_basis_matrix_B, eigen_values, error_estimate] = ...
    pod_extension( snapshot_matrix, params, min_eigen_value, ...
    max_reduced_basis, inner_product_matrix )
%POD_EXTENSION Summary of this function goes here
%   Detailed explanation goes here

ns = size(snapshot_matrix,2);
M = inner_product_matrix;

[ V, D] = eig(snapshot_matrix' * M * snapshot_matrix);
[d,ind] = sort(diag(D),'descend');
Ds = D(ind,ind);
eigen_values = d;
Vs = V(:,ind);

if max_reduced_basis > sum(d >= min_eigen_value)
    max_reduced_basis = sum(d >= min_eigen_value);
    disp('Size of reduced basis space revised based on minimum eigenvalue criteria');
end

Ds = Ds(1:max_reduced_basis,1:max_reduced_basis);
Vs = Vs(:,1:max_reduced_basis);
%R = [eye(max_reduced_basis);...
%    zeros(ns - max_reduced_basis,max_reduced_basis)];

reduced_basis_matrix_B = zeros(params.ndofs,max_reduced_basis);
for i = 1:1:max_reduced_basis
    reduced_basis_matrix_B(:,i) = snapshot_matrix * Vs(:,i) ...
        / sqrt((snapshot_matrix * Vs(:,i))' * M * ...
        (snapshot_matrix * Vs(:,i)));
end

error_estimate = sqrt(sum(eigen_values(...
    size(reduced_basis_matrix_B,2)+1:size(snapshot_matrix,2)))) / sqrt(sum(eigen_values));

end