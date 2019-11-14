function [ res, res1, res2, res3 ] = assemble_bilinear_side_extension...
    ( params, paramsP, grid, mu, c11 )
%ASSEMBLE_BILINEAR_SIDE_EXTENSION Summary of this function goes here
%   Detailed explanation goes here

diffusive_term = ldg_evaluate_basis_derivative_assembly_extension...
    ( params, grid);

% % TODO : Using old assembly. Use extension file after correction.
% coercivity_term = n_u_tensor_jump_n_phi_tensor_jump_assembly...
%     (grid,params);
coercivity_term = n_u_tensor_jump_n_phi_tensor_jump_assembly_offline...
    ( params, grid);

flux_approximation = del_phi_average_n_phi_tensor_jump_assembly_offline...
    ( params, grid );

close all

res1 = diffusive_term;
res2 = coercivity_term;
res3 = flux_approximation;

res = mu * (diffusive_term{1} + diffusive_term{4}) + ...
    c11 * coercivity_term.res - ...
    mu * (flux_approximation.res) - mu * (flux_approximation.res)';

end