function [ res ] = n_u_tensor_jump_n_phi_tensor_jump_plus_minus_integral_12...
    ( params, grid, tria_index, local_vertex_index, qdeg )
%N_U_TENSOR_JUMP_N_PHI_TENSOR_JUMP_PLUS_PLUS_INTEGRAL Summary of this function goes here
%   Detailed explanation goes here

if nargin == 4
    qdeg = params.qdeg;
end

f = @(llcoord) n_u_tensor_jump_n_phi_tensor_jump_plus_minus_local_12...
    ( llcoord, params, grid, tria_index, local_vertex_index);

res = intervalquadrature(qdeg,f) * grid.EL(tria_index,local_vertex_index);

end