function [ res ] = del_u_average_n_v_tensor_jump_plus_minus_local( llcoord,...
    params,grid,tria_index,local_vertex_index)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

res = zeros(params.ndofs_per_element);
JIT = [grid.JIT(tria_index,:,1)',grid.JIT(tria_index,:,2)']; 
tria_index_neighbour = grid.NBI(tria_index,local_vertex_index);
if tria_index_neighbour > 0
    local_vertex_index_neighbour = ...
        find(grid.NBI(tria_index_neighbour,:) == tria_index);
end
N = [grid.NX(tria_index,local_vertex_index) ...
    grid.NY(tria_index,local_vertex_index)];
N_neighbour = -N;
velocity_basis_gradient_self = ldg_evaluate_basis_derivative(...
    llocal2local(grid,local_vertex_index,llcoord),params);
velocity_basis_gradient_self = cell2mat(velocity_basis_gradient_self);
velocity_basis_gradient_self = sum(velocity_basis_gradient_self,1);
velocity_basis_gradient_self = reshape(velocity_basis_gradient_self,...
    1,2,params.ndofs_per_element);
%2 is because of 2-d grid
% now velocity_basis_gradient_self(:,:,j) represents gradient of jth theoretical
% basis function
velocity_basis_neighbour = ldg_evaluate_basis...
    (llocal2local(grid,local_vertex_index_neighbour,llcoord),params);
velocity_basis_neighbour = sum(velocity_basis_neighbour,2);
%velocity_basis_neighbour = reshape(velocity_basis_neighbour,...
%    params.dimrange,params.ndofs_per_element/params.dimrange);
%velocity_basis_neighbour = velocity_basis_neighbour';
%now each row represents one theoretical basis function
if tria_index_neighbour > 0
    for a=1:1:params.dimrange
        for i = a:params.dimrange:params.ndofs_per_element
            for j = a:params.dimrange:params.ndofs_per_element
                res(i,j) = res(i,j) +...
                    0.5 * (velocity_basis_neighbour(i) * N_neighbour) * ...
                    (velocity_basis_gradient_self(:,:,j)*JIT')';
            end
        end
    end
else
    disp('Neighour evaluation for dirichlet boundary, returning zeros')
end
end
%check formula