function [ res ] = ldg_evaluate_basis_derivative_assembly_extension_online...
    ( offline_assembly, para_mapping, params, grid, el_subd)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

gids = ldg_global_dof_index(params,grid);

res = zeros(params.ndofs);

for i=1:1:length(el_subd)
    determinant = 1 / det(para_mapping{i});
    affine_multiplier_11 = para_mapping{i}(1,:) * ...
        para_mapping{i}(1,:)';
    affine_multiplier_12 = para_mapping{i}(1,:) * ...
        para_mapping{i}(2,:)';
    affine_multiplier_21 = para_mapping{i}(2,:) * ...
        para_mapping{i}(1,:)';
    affine_multiplier_22 = para_mapping{i}(2,:) * ...
        para_mapping{i}(2,:)';
    ids = gids(el_subd{i},:);
    res(ids,ids) = determinant * (res(ids,ids) + ...
        affine_multiplier_11 * offline_assembly{1}(ids,ids) + ...
        affine_multiplier_12 * offline_assembly{2}(ids,ids) + ...
        affine_multiplier_21 * offline_assembly{3}(ids,ids) + ...
        affine_multiplier_22 * offline_assembly{4}(ids,ids));
end

if params.show_sparsity == true
    figure()
    spy(res) % visualise sparsity pattern
    title('Spy of \nabla \phi : \nabla \phi')
    axis equal
    axis tight
    disp('Observe all graphs of \nabla \phi : \nabla \phi')
    pause()
end

end