function [ res ] = assemble_lhs_continuity( params, paramsP, grid, qdeg )
%ASSEMBLE_LHS_CONTINUITY Summary of this function goes here
%   Detailed explanation goes here

res1 = pressure_velocity_continuity_assembly( params, paramsP, grid, qdeg);
res1 = sparse(res1);
res2 = pressure_average_velocity_basis_jump_assembly(grid,params,paramsP);
res2 = sparse(res2.res);

res = res1 + res2;

end

