function [ res ] = n_u_tensor_jump_n_phi_tensor_jump_minus_plus_local_12...
    ( llcoord, params, grid, tria_index, local_vertex_index)
%N_U_TENSOR_JUMP_N_PHI_TENSOR_JUMP_PLUS_PLUS_LOCAL Summary of this function goes here
%   Detailed explanation goes here

res = zeros(params.ndofs_per_element);
tria_index_neighbour = grid.NBI(tria_index, local_vertex_index);
N = [grid.NX(tria_index,local_vertex_index) ...
    grid.NY(tria_index,local_vertex_index)];
N_neighbour = -N;
velocity_basis_self = ldg_evaluate_basis...
    (llocal2local(grid,local_vertex_index,llcoord), params);

if tria_index_neighbour > 0
    local_vertex_index_neighbour = ...
        find(grid.NBI(tria_index_neighbour,:)==tria_index);
    velocity_basis_neighbour = ldg_evaluate_basis(llocal2local(grid,...
        local_vertex_index_neighbour,llcoord), params);
    for a = 1:1:params.dimrange
        for i = a:params.dimrange:params.ndofs_per_element
            for j = a:params.dimrange:params.ndofs_per_element
                res(i,j) = res(i,j) + ...
                    (N_neighbour(1) * sum(velocity_basis_neighbour(j,:))) * ...
                    (N(2) * sum(velocity_basis_self(i,:)));
            end
        end
    end
else
    % Do nothing as \phi^- is zero.
end

end