function [ res ] = ldg_evaluate_basis_derivative_integral...
    ( params, k, grid, qdeg )
%LDG_EVALUATE_BASIS_DERIVATIVE_LOCAL_INTEGRAL Summary of this function goes here
%   Detailed explanation goes here

if nargin == 3
    qdeg = params.qdeg;
end

f=@(lcoord) ldg_evaluate_basis_derivative_local(lcoord, params, k, grid);
res=triaquadrature(qdeg,f)*(2*grid.A(k));%2*grid.A(k) is jacobian of element k;

%spy(res); % Visualisation of sparsity pattern

end

