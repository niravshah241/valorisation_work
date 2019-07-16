function [ res ] = source_assembly_extension_online...
    ( offline_vector, params, grid, para_mapping, el_subd)
%SOURCE_ASSEMBLY_EXTENSION_ONLINE Summary of this function goes here
%   Detailed explanation goes here

gids = ldg_global_dof_index(params,grid);
res = zeros(params.ndofs,1);

for i = 1:1:length(el_subd)
    res(gids(el_subd{i},:),1) = offline_vector(gids(el_subd{i},:),1) / ...
        det(para_mapping{i});
end

end