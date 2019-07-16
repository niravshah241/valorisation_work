function [ res, res1, res2, res3, res4 ] = assemble_linear_side...
    ( grid, params, paramsP, qdeg, mu, c11)
%ASSEMBLE_LINEAR_SIDE Summary of this function goes here
%   Detailed explanation goes here


res1 = source_assembly(params, paramsP, grid);
res1 = sparse(res1);
res2 = t_v_assembly( params, paramsP, grid, qdeg );
res2 = sparse(res2);
res3 = u_d_v_assembly( grid, params, qdeg);
res3 = sparse(res3);
res4 = n_uD_tensor_del_v_assembly( params, grid, qdeg);
res4 = sparse(res4);

rhs = res1 + res2 + c11 * res3 - mu * res4;

res = sparse(rhs);

close all

end