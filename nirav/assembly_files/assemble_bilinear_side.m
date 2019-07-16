function [ res, res1, res2, res3 ] = assemble_bilinear_side...
    ( grid, params, paramsP, qdeg, mu, c11)
%ASSEMBLE_BILINEAR_SIDE Summary of this function goes here
%   Detailed explanation goes here


res1 = ldg_evaluate_basis_derivative_assembly(params, grid );
res1 = sparse(res1);
res2 = n_u_tensor_jump_n_phi_tensor_jump_assembly(grid,params);
res2 = sparse(res2.res);
res3 = del_u_average_n_v_tensor_jump_assembly(grid,params);
res3 = sparse(res3.res);

res = mu * res1 + c11 * res2 - mu * res3 - mu * res3';

end