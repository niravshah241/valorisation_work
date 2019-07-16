function [ params, paramsP, rhs, stifness_matrix] = ...
    assemble_stifness_matrix...
    ( params, paramsP, grid, qdeg, mu, c11 )
%variational form assembly
% Momentum : AU + BP = C
% Continuity      B'U = D

[params.bilinear_side, params.bilinear_res1, params.bilinear_res2, ...
    params.bilinear_res3] = assemble_bilinear_side(grid, params,...
    paramsP, qdeg, mu, c11); %A
params.bilinear_side_pressure_terms =...
    assemble_bilinear_side_pressure_terms( params, paramsP, grid);%B
[params.linear_side, params.linear_res1, params.linear_res2, ...
params.linear_res3, params.linear_res4] = ...
assemble_linear_side(grid, params, paramsP, qdeg, mu, c11);%C

params.rhs_continuity = assemble_rhs_continuity( params, paramsP, grid, qdeg);%D
params.lhs_continuity = assemble_lhs_continuity( params, paramsP, grid, qdeg);%B'

dofs = [params.dofs;paramsP.dofs];
linear_stifness_matrix = params.linear_side; 
continuity_pressure = zeros(paramsP.ndofs); 
bilinear_stifness_matrix = [params.bilinear_side,params.bilinear_side_pressure_terms];
continuity_stifness_matrix = [params.lhs_continuity,continuity_pressure];
stifness_matrix = [bilinear_stifness_matrix;continuity_stifness_matrix];
rhs = [params.linear_side;params.rhs_continuity];

%converting to sparse matrix

stifness_matrix = sparse(stifness_matrix);

params.bilinear_res1 = sparse(params.bilinear_res1);
params.bilinear_res2 = sparse(params.bilinear_res2);
params.bilinear_res3 = sparse(params.bilinear_res3);

params.linear_res1 = sparse(params.linear_res1);
params.linear_res2 = sparse(params.linear_res2);
params.linear_res3 = sparse(params.linear_res3);
params.linear_res4 = sparse(params.linear_res4);

end