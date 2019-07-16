function [ res ] = ldg_evaluate_basis_derivative_integral_11_extension...
    ( params, k, grid, qdeg )

if nargin == 3
    qdeg = params.qdeg;
end

f = @(lcoord) ldg_evaluate_basis_derivative_local_11_extension...
    (lcoord, params, k, grid);
res = triaquadrature(qdeg,f) * (2 * grid.A(k));
%2*grid.A(k) is jacobian of element k;

%spy(res); % Visualisation of sparsity pattern

end

