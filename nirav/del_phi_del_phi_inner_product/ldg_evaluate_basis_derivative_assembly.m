function res = ldg_evaluate_basis_derivative_assembly(params, grid, qdeg)

if nargin == 2
    qdeg = params.qdeg;
end

% getting global ids of dof
gids = ldg_global_dof_index(params, grid);

% Allocating space for output
res=zeros(params.ndofs,params.ndofs);

%Evaluation and assembly
for k=1:1:grid.nelements
    ids=gids(k,:); % Extracting global dof ids of current element
    % tmp is temporary variable to store current values out of integral,
    % and is overwritten in each loop
    tmp = ldg_evaluate_basis_derivative_integral( params, k, grid);
    res(ids,ids)=res(ids,ids)+tmp;
end

if params.show_sparsity == true
    figure()
    spy(res) % visualise sparsity pattern
    title('Spy of \nabla \phi : \nabla \phi')
    axis equal
    axis tight
    disp('Observe all graphs of \nabla \phi : \nabla \phi')
    pause()
end

res = (res+res')/2;

end