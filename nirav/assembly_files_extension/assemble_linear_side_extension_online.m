function [ res ] = assemble_linear_side_extension_online...
    ( params, paramsP, grid, source_vector_offline, ...
    offline_rhs_vector, el_subd, para_mapping )
%ASSEMBLE_LINEAR_SIDE_EXTENSION_ONLINE Summary of this function goes here
%   Detailed explanation goes here

source_vector_new = source_assembly_extension_online...
    ( source_vector_offline, params, grid, para_mapping, el_subd);

res = offline_rhs_vector - source_vector_offline + source_vector_new;
res = sparse(res);

end

