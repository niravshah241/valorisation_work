function [ res ] = pressure_velocity_continuity_local_12_extension...
    ( lcoord, params, paramsP, k, grid)
%PRESSURE_VELOCITY_CONTINUITY_LOCAL_11_EXTENSION Summary of this function goes here
%   Detailed explanation goes here

JIT=[grid.JIT(k,:,1)',grid.JIT(k,:,2)'];

pressure_basis = ldg_evaluate_basis(lcoord,paramsP);

basis_derivative = ldg_evaluate_basis_derivative(lcoord,params);

JIT = [grid.JIT(k,:,1)',grid.JIT(k,:,2)']; % Jacobian Inverse Transpose

for k = 1:1:params.ndofs_per_element
    velocity_basis_derivative{k} = (JIT * (basis_derivative{k})')';
end

%velocity_basis_derivative has basis function along row i.e. same as
%local_basis_derivative

res = zeros(paramsP.ndofs_per_element,params.ndofs_per_element);

for i = 1:1:paramsP.ndofs_per_element
    for j = 1:1:params.ndofs_per_element
        res(i,j) = res(i,j) + pressure_basis(i,:) * ...
            velocity_basis_derivative{j}(1,2);
    end
end

end