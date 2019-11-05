function [ reduced_basis_matrix_B, eigen_values, error_estimate] = ...
    pod_extension_velocity( snapshot_matrix, params, min_eigen_value, ...
    max_reduced_basis, inner_product_matrix )
%POD_EXTENSION Summary of this function goes here
%   Detailed explanation goes here

error_estimate = zeros(params.dimrange,1);

for k = 1:1:params.dimrange
    ns = size(snapshot_matrix,2);
    M = inner_product_matrix...
        (k:params.dimrange:params.ndofs,k:params.dimrange:params.ndofs);
    [ V, D] = eig((snapshot_matrix(k:params.dimrange:params.ndofs,:))' * ...
        M * snapshot_matrix(k:params.dimrange:params.ndofs,:));
    [d,ind] = sort(diag(D),'descend');
    Ds = D(ind,ind);
    eigen_values{k} = d;
    Vs = V(:,ind);
    
    if max_reduced_basis(k) > sum(d >= min_eigen_value(k))
        max_reduced_basis(k) = sum(d >= min_eigen_value(k));
        disp('Size of reduced basis space revised based on minimum eigenvalue criteria');
    end
    
    Ds = Ds(1:max_reduced_basis(k),1:max_reduced_basis(k));
    Vs = Vs(:,1:max_reduced_basis(k));
    %R = [eye(max_reduced_basis);...
    %    zeros(ns - max_reduced_basis,max_reduced_basis)];
    
    reduced_basis_matrix_B{k} = zeros(params.ndofs/params.dimrange,...
        max_reduced_basis(k));
    for i = 1:1:max_reduced_basis(k)
        reduced_basis_matrix_B{k}(:,i) = ...
            snapshot_matrix(k:params.dimrange:params.ndofs,:) * Vs(:,i) ...
            / sqrt((snapshot_matrix(k:params.dimrange:params.ndofs,:) * ...
            Vs(:,i))' * M * ...
            (snapshot_matrix (k:params.dimrange:params.ndofs,:) * Vs(:,i)));
    end
    
    error_estimate(k) = sqrt(sum(eigen_values{k}(...
        size(reduced_basis_matrix_B{k},2)+1:size(snapshot_matrix,2)))) / ...
        sqrt(sum(eigen_values{k}));
end

end