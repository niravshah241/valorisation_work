function [ res ] = pressure_velocity_continuity_assembly_extension_online...
    ( offline_assembly, para_mapping, params, paramsP, grid, el_subd)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

gids_v = ldg_global_dof_index(params,grid);
gids_p = ldg_global_dof_index(paramsP,grid);

res = zeros( paramsP.ndofs, params.ndofs);

for i=1:1:length(el_subd)
    determinant = 1 / det(para_mapping{i});
    affine_multiplier_11 = para_mapping{i}(1,1);
    affine_multiplier_12 = para_mapping{i}(2,1);
    affine_multiplier_21 = para_mapping{i}(1,2);
    affine_multiplier_22 = para_mapping{i}(2,2);
    
%     %%Check whether below is correct instead of above
%     affine_multiplier_11 = para_mapping{i}(1,1);
%     affine_multiplier_12 = para_mapping{i}(2,1);
%     affine_multiplier_21 = para_mapping{i}(1,2);
%     affine_multiplier_22 = para_mapping{i}(2,2);
    
    ids_v = gids_v(el_subd{i},:);
    ids_p = gids_p(el_subd{i},:);
    
    res(ids_p,ids_v) = determinant * (res(ids_p,ids_v) + ...
        affine_multiplier_11 * offline_assembly{1}(ids_p,ids_v) + ...
        affine_multiplier_12 * offline_assembly{2}(ids_p,ids_v) + ...
        affine_multiplier_21 * offline_assembly{3}(ids_p,ids_v) + ...
        affine_multiplier_22 * offline_assembly{4}(ids_p,ids_v));
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