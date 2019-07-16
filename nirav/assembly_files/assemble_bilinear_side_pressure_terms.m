function [ res ] = assemble_bilinear_side_pressure_terms...
    ( params, paramsP, grid)
%ASSEMBLE_BILINEAR_SIDE_PRESSURE_TERMS Summary of this function goes here
%   Detailed explanation goes here

res1 = pressure_velocity_continuity_assembly( params, paramsP, grid);
res1 = sparse(res1);
res2 = pressure_average_velocity_basis_jump_assembly(grid,params,paramsP);
res2 = sparse(res2.res);


res = res1' + res2';

end

