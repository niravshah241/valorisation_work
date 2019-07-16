function [ res ] = n_u_tensor_jump_n_phi_tensor_jump_plus_minus_local( llcoord,...
    params, grid, tria_index, local_vertex_index)
%N_U_TENSOR_JUMP_N_PHI_TENSOR_JUMP_PLUS_PLUS_LOCAL Summary of this function goes here
%   Detailed explanation goes here

res = zeros(params.ndofs_per_element);
tria_index_neighbour = grid.NBI(tria_index, local_vertex_index);
N = [grid.NX(tria_index,local_vertex_index) grid.NY(tria_index,local_vertex_index)];
N_neighbour = -N;
velocity_basis_self = ldg_evaluate_basis(llocal2local(grid,local_vertex_index,...
    llcoord), params);
velocity_basis_self = sum(velocity_basis_self,2);
%velocity_basis_self_theory = reshape(velocity_basis_self_theory,params.dimrange,...
%  params.ndofs_per_element/params.dimrange);
%velocity_basis_self_theory = velocity_basis_self_theory';
%now every row corressponds to one basis function
if tria_index_neighbour > 0
    local_vertex_index_neighbour = find(grid.NBI(tria_index_neighbour,:)==tria_index);
    velocity_basis_neighbour = ldg_evaluate_basis(llocal2local(grid,...
        local_vertex_index_neighbour,llcoord), params);
    velocity_basis_neighbour = sum(velocity_basis_neighbour,2);
    %velocity_basis_neighbour_theory = reshape(velocity_basis_neighbour_theory,...
    %   params.dimrange,params.ndofs_per_element/params.dimrange);
    %velocity_basis_neighbour_theory = velocity_basis_neighbour_theory';
    % now every row corressponds to one basis function
    for a = 1:1:params.dimrange
        for i = a:params.dimrange:params.ndofs_per_element
            for j = a:params.dimrange:params.ndofs_per_element
                res(i,j) = res(i,j) + ...
                    (N*velocity_basis_self(j))*...
                    (N_neighbour*velocity_basis_neighbour(i))';
                % 2 is for 2-d grid
            end
        end
    end
else
    disp('Dirichlet boundary for neighbour part, returning zeros')
end


end