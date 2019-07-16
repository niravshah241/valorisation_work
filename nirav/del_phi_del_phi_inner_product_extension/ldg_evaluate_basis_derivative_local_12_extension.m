function res = ldg_evaluate_basis_derivative_local_12_extension...
    ( lcoord, params, k, grid)
%LDG_EVALUATE_BASIS_DERIVATIVE_LOCAL: assembly of \phi:\phi where \phi is
%global basis function for velocity. \nabla \phi = JIT * \hat{\nabla \phi} where
%\hat{\phi} is local basis function for velocity
%Params needs fields: dimrange, pdeg, ndofs, ndofs_per_element
%k is element number and grid should have field JIT (Jacobian Inverse Transposed)

basis_derivative = ldg_evaluate_basis_derivative(lcoord,params);

JIT = [grid.JIT(k,:,1)',grid.JIT(k,:,2)']; % Jacobian Inverse Transpose

for k = 1:1:params.ndofs_per_element
    global_basis_derivative{k} = (JIT * (basis_derivative{k})')';
end

%global_basis_derivative has basis function along row i.e. same as
%local_basis_derivative

res = zeros(params.ndofs_per_element, params.ndofs_per_element);

for k = 1:1:params.dimrange
    for i = k:params.dimrange:params.ndofs_per_element
        for j = k:params.dimrange:params.ndofs_per_element
            %res(i,j) = sum(global_basis_derivative{i},1) *...
             %   params.mu * sum(global_basis_derivative{j},1)';
             basis_i = sum(global_basis_derivative{i},1);
             basis_j = sum(global_basis_derivative{j},1);
             res(i,j) = basis_i(1) * basis_j(2);
        end
    end
end

end