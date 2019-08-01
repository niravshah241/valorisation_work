function [ params_reduced, paramsP_reduced, reduced_rhs, ...
    reduced_stiffness_matrix] = galerkin_extension( params, paramsP, ...
    grid, rhs, stiffness_matrix, reduced_basis_matrix_B_pressure, ...
    reduced_basis_matrix_B_velocity);
%GALERKIN_EXTENSION Summary of this function goes here
%   Detailed explanation goes here
params_reduced = params;
paramsP_reduced = paramsP;

rhs_reduced = zeros(red_dim_velocity + red_dim_pressure,1);
rhs_reduced(1:red_dim_velocity) = B_velocity' * rhs(1:params.ndofs);
rhs_reduced(red_dim_velocity + 1:red_dim_velocity + red_dim_pressure) = ...
    B_pressure' * rhs(params.ndofs + 1,params.ndofs + paramsP.ndofs);
stifness_matrix_reduced = zeros(red_dim_velocity+red_dim_pressure);
stifness_matrix_reduced(1:red_dim_velocity,1:red_dim_velocity) = ...
    B_velocity' * stifness_matrix(1:params.ndofs,1:params.ndofs) * ...
    B_velocity;
stifness_matrix_reduced(1:red_dim_velocity,red_dim_velocity+1:...
    red_dim_velocity + red_dim_pressure) = B_velocity' * ...
    stifness_matrix(1:params.ndofs,params.ndofs + 1 : ...
    params.ndofs + paramsP.ndofs) * B_pressure;
stifness_matrix_reduced(red_dim_velocity + 1:red_dim_velocity + ...
    red_dim_pressure, 1:red_dim_velocity) = B_pressure' * ...
    stifness_matrix(params.ndofs + 1 : params.ndofs + paramsP.ndofs, ...
    1:params.ndofs) * B_velocity;
stifness_matrix_reduced = sparse(stifness_matrix_reduced);

[ params_reduced, paramsP_reduced] = solve_plot_solution_schur...
    ( params_reduced, paramsP_reduced, grid, rhs_reduced, ...
    stifness_matrix_reduced, 0);
params_reduced.dofs = reduced_basis_matrix_B_velocity * ...
    params_reduced.dofs;
paramsP_reduced.dofs = reduced_basis_matrix_B_pressure * ...
    paramsP_reduced.dofs;

for i=1:1:params.dimrange
    figure()
    axis equal
    [scalar_dofs, scalar_df_info] = ...
        ldg_scalar_component(params.dofs,i,params);
    sdf = ldgdiscfunc(scalar_dofs,scalar_df_info);
    disp(['Plotting ',num2str(i),' degree of freedom (Schur)'])
    %subplot(params.dimrange,1,i)
    %title(['Velocity degree of freedom number ',num2str(i)])
    if i==1
        title(['Velocity in x direction (Schur)'])
    else
        title(['Velocity in y direction (Schur)'])
    end
    axis equal
    axis tight
    ldg_plot(sdf,grid,params);
    plot(grid);
end

for i=1:1:paramsP.dimrange
    figure()
    [scalar_dofs, scalar_df_info] = ...
        ldg_scalar_component(paramsP.dofs,i,paramsP);
    sdf = ldgdiscfunc(scalar_dofs,scalar_df_info);
    disp(['Plotting ',num2str(i),' degree of freedom (for pressure)'])
    %subplot(paramsP.dimrange,1,i)
    title('Pressure (Schur)')
    axis equal
    axis tight
    ldg_plot(sdf,grid,paramsP);
    plot(grid);
end

end