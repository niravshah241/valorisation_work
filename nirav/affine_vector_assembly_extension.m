function [ affine_vector_assembly ] = ...
    affine_vector_assembly_extension...
    ( vector, affine_factors, el_subd, params, grid )
%AFFINE_ASSEMBLY_EXTENSION Summary of this function goes here
%   Detailed explanation goes here

gids = ldg_global_dof_index(params,grid);
affine_vector_assembly = ones(length(vector),1);

for k = 1:1:length(el_subd)
    elements = el_subd{k};
    for i = 1:1:length(elements)
        affine_vector_assembly(gids(elements(i),:)) = ...
            affine_factors(k) * vector(gids(elements(i),:));
    end
end


% for i = 1:1:grid.nelements
%     k = get_subdomain_number_extension(i,el_subd);
%     affine_vector_assembly(gids(i,:)) = ...
%         affine_factors(k) * vector(gids(i,:));
% end
end