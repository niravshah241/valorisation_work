function [ res ] = del_u_average_n_v_tensor_jump_plus_plus_local_21...
    ( llcoord, params, grid, tria_index, local_vertex_index)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

res = zeros(params.ndofs_per_element);
JIT = [grid.JIT(tria_index,:,1)',grid.JIT(tria_index,:,2)'];
tria_index_neighbour = grid.NBI(tria_index,local_vertex_index);
N = [grid.NX(tria_index,local_vertex_index) ...
    grid.NY(tria_index,local_vertex_index)];

basis_derivative = ldg_evaluate_basis_derivative(...
    llocal2local(grid,local_vertex_index,llcoord),params);

velocity_basis_self = ldg_evaluate_basis...
    (llocal2local(grid,local_vertex_index,llcoord),params);

for k = 1:1:params.ndofs_per_element
    global_basis_derivative{k} = (JIT * (basis_derivative{k})')';
end

if tria_index_neighbour > 0
    for a=1:1:params.dimrange
        for i = a:params.dimrange:params.ndofs_per_element
            for j = a:params.dimrange:params.ndofs_per_element
                temp = sum(global_basis_derivative{j},1);
                res(i,j) = res(i,j) + 0.5 * sum(velocity_basis_self(i,:)) * ...
                    N(2) * temp(1);
            end
        end
    end
else
    for a=1:1:params.dimrange
        for i = a:params.dimrange:params.ndofs_per_element
            for j = a:params.dimrange:params.ndofs_per_element
                temp = sum(global_basis_derivative{j},1);
                res(i,j) = res(i,j) + sum(velocity_basis_self(i,:)) * ...
                    N(2) * temp(1);
            end
        end
    end
end
end