function [ rhs, params, linear_side, source_vector, ...
    linear_side_continuity] = assemble_rhs_extension...
    ( params, paramsP, grid, mu, c11, qdeg)
%ASSEMBLE_RHS_EXTENSION Summary of this function goes here
%   Detailed explanation goes here

if nargin == 5
    qdeg = params.qdeg;
end

rhs = zeros(params.ndofs + paramsP.ndofs,1);
[ linear_side, source_vector ] = ...
    assemble_linear_side_extension( params, paramsP, grid, mu, c11 );
linear_side_continuity = assemble_rhs_continuity( params, paramsP, grid);
rhs(1:params.ndofs) = linear_side;
rhs(params.ndofs+1:params.ndofs+paramsP.ndofs) = linear_side_continuity;
params.linear_side = rhs(1:params.ndofs);
params.rhs_continuity = rhs(params.ndofs+1:params.ndofs+paramsP.ndofs);

end