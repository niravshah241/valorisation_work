function [ res ] = del_u_average_n_v_tensor_jump_minus_plus_local( llcoord,...
    params,grid,tria_index,local_vertex_index)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

res = zeros(params.ndofs_per_element);
tria_index_neighbour = grid.NBI(tria_index,local_vertex_index);
JIT = [grid.JIT(tria_index_neighbour,:,1)' ...
    grid.JIT(tria_index_neighbour,:,2)']; 
if tria_index_neighbour > 0
    local_vertex_index_neighbour = ...
        find(grid.NBI(tria_index_neighbour,:) == tria_index);
end
N = [grid.NX(tria_index,local_vertex_index) ...
    grid.NY(tria_index,local_vertex_index)];
N_neighbour = -N;
velocity_basis_gradient_neighbour = ldg_evaluate_basis_derivative(...
    llocal2local(grid,local_vertex_index_neighbour,llcoord),params);
velocity_basis_gradient_neighbour = cell2mat(velocity_basis_gradient_neighbour);
velocity_basis_gradient_neighbour = sum(velocity_basis_gradient_neighbour,1);
velocity_basis_gradient_neighbour = reshape(velocity_basis_gradient_neighbour,...
    1,2,params.ndofs_per_element);
%2 is because of 2-d grid
% now velocity_basis_gradient_self(:,:,j) represents gradient of jth theoretical
% basis function
velocity_basis_self = ldg_evaluate_basis...
    (llocal2local(grid,local_vertex_index,llcoord),params);
velocity_basis_self = sum(velocity_basis_self,2);
%velocity_basis_self = reshape(velocity_basis_self,...
%    params.dimrange,params.ndofs_per_element/params.dimrange);
%velocity_basis_self = velocity_basis_self';
%now each row represents one theoretical basis function
if tria_index_neighbour > 0
    for a = 1:1:params.dimrange
        for i = a:params.dimrange:params.ndofs_per_element
            for j = a:params.dimrange:params.ndofs_per_element
                res(i,j) = res(i,j) +...
                    0.5 * (velocity_basis_self(i) * N) * ...
                    (velocity_basis_gradient_neighbour(:,:,j) * JIT')';
            end
        end
    end
else
    disp('Neighour evaluation for dirichlet boundary, returnin zeros')
end
end
%check formula