function [ stiffness_matrix_online, params ] = ...
    assemble_stiffness_matrix_online_extension( params, paramsP, ...
    grid, diffusive_term, coercivity_term, ...
    flux_approximation, para_mapping, el_subd, transformed_grid, c11, ...
    mu, incompressibility_term, pressure_flux_term)
%ASSEMBLE_STIFFNESS_MATRIX_ONLINE_EXTENSION Summary of this function goes here
%   Detailed explanation goes here

[ bilinear_side_online ] = assemble_bilinear_side_extension_online...
    ( params, paramsP, grid, diffusive_term, coercivity_term, ...
    flux_approximation, para_mapping, el_subd, transformed_grid, c11, mu);
[ bilinear_side_pressure_terms_online ] = ...
    assemble_bilinear_side_pressure_terms_extension_online...
    ( params, paramsP, grid, incompressibility_term, ...
    pressure_flux_term, para_mapping, el_subd);
stiffness_matrix_online = sparse(zeros(params.ndofs+paramsP.ndofs));
stiffness_matrix_online(1:params.ndofs,1:params.ndofs) = ...
    bilinear_side_online;
stiffness_matrix_online(1:params.ndofs,...
    params.ndofs+1:params.ndofs+paramsP.ndofs) = ...
    bilinear_side_pressure_terms_online;
stiffness_matrix_online(params.ndofs+1:params.ndofs+paramsP.ndofs,...
    1:params.ndofs) = bilinear_side_pressure_terms_online';
params.bilinear_side = bilinear_side_online;
params.bilinear_side_pressure_terms = ...
    bilinear_side_pressure_terms_online;

end