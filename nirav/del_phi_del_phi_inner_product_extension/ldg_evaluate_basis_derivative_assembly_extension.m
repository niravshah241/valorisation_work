function res = ldg_evaluate_basis_derivative_assembly_extension...
    ( params, grid)

if nargin == 2
    qdeg = params.qdeg;
end

gids = ldg_global_dof_index(params,grid);
%res = zeros(params.ndofs,params.ndofs);
res_11 = sparse(zeros(params.ndofs,params.ndofs));
res_12 = sparse(zeros(params.ndofs,params.ndofs));
res_21 = sparse(zeros(params.ndofs,params.ndofs));
res_22 = sparse(zeros(params.ndofs,params.ndofs));

for k = 1:1:grid.nelements
    tmp_11 = ldg_evaluate_basis_derivative_integral_11_extension...
        (params, k, grid);
    tmp_12 = ldg_evaluate_basis_derivative_integral_12_extension...
        (params, k, grid);
    tmp_21 = ldg_evaluate_basis_derivative_integral_21_extension...
        (params, k, grid);
    tmp_22 = ldg_evaluate_basis_derivative_integral_22_extension...
        (params, k, grid);
    ids = gids(k,:);
    %res(ids,ids) = res(ids,ids) + tmp;
    res_11(ids,ids) = res_11(ids,ids) + tmp_11;
    res_12(ids,ids) = res_12(ids,ids) + tmp_12;
    res_21(ids,ids) = res_21(ids,ids) + tmp_21;
    res_22(ids,ids) = res_22(ids,ids) + tmp_22;
end

res{1} = res_11;
res{2} = res_12;
res{3} = res_21;
res{4} = res_22;

%res = (res+res')/2;

end