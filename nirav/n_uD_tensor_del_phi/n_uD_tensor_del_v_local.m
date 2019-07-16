function [ res ] = n_uD_tensor_del_v_local( llcoord, params, grid,...
    tria_index, local_vertex_index )
%N_UD_TENSOR_DEL_V_LOCAL Summary of this function goes here
%   Detailed explanation goes here

JIT = [grid.JIT(tria_index,:,1)',grid.JIT(tria_index,:,2)']; 

res = zeros(params.ndofs_per_element,1);
lcoord = llocal2local(grid,local_vertex_index,llcoord);
glob = local2global(grid,tria_index,lcoord,params);
N = [grid.NX(tria_index,local_vertex_index) grid.NY(tria_index,local_vertex_index)];
uD = dirichlet_values(glob,params);
uD_N = kron(uD,N);
uD_N = reshape(uD_N,1,4);
velocity_basis_gradient = ldg_evaluate_basis_derivative(lcoord,params);
velocity_basis_gradient = cell2mat(velocity_basis_gradient);
velocity_basis_gradient_reshaped =...
    reshape(velocity_basis_gradient,params.dimrange,2,params.ndofs_per_element);
%2 is for triagrid
for i=1:1:params.ndofs_per_element
    velocity_basis_global = velocity_basis_gradient_reshaped(:,:,i)*JIT';
    temp = reshape(velocity_basis_global,1,4);
    res(i) = res(i) + temp * uD_N';
end

% res = zeros(params.ndofs_per_element,1);
% 
% lcoord = llocal2local(grid,local_vertex_index,llcoord);
% glob = local2global(grid,tria_index,lcoord,params);
% uD = dirichlet_values(glob);
% N = [grid.NX(tria_index,local_vertex_index)...
%     grid.NY(tria_index,local_vertex_index)];
% velocity_basis = ldg_evaluate_basis(lcoord,params);
% velocity_basis_gradient = ldg_evaluate_basis_derivative(lcoord,params);
% velocity_basis_gradient = cell2mat(velocity_basis_gradient);
% velocity_basis_gradient = sum(velocity_basis_gradient,1);
% velocity_basis_gradient = reshape(velocity_basis_gradient,...
%     1,2,params.ndofs_per_element);
% %2 is because of 2-d grid
% % now velocity_basis_gradient(:,:,j) represents gradient of jth theoretical
% % basis function
% for i=1:1:params.ndofs_per_element
%     %a = rem(i,params.dimrange);
%     %if a == 0
%     %    a=params.dimrange;
%     %end
%     res(i) = res(i) + (uD(i)*N)*(velocity_basis_gradient(:,:,i))';
% end

end