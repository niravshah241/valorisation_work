tic();
[ transformed_grid, F_transformation_matrix, ...
    C_translation_vector] = transform_grid( params, grid, mu_x, mu_y);
close all

for i = 1:1:length(F_transformation_matrix)
    para_mapping{i} = inv(F_transformation_matrix{i});
end

assert(length(para_mapping) == length(el_subd),...
    'Number of subdomains and number of parametric mappings not same');

[ stiffness_matrix_online, params ] = ...
    assemble_stiffness_matrix_online_extension( params, paramsP, ...
    grid, diffusive_term, coercivity_term, ...
    flux_approximation, para_mapping, el_subd, transformed_grid, c11, ...
    mu, incompressibility_term, pressure_flux_term);

[ rhs_online, params ] = assemble_rhs_online_extension...
    ( params, paramsP, grid, source_vector_offline, ...
    linear_side_offline, el_subd, para_mapping, linear_side_continuity);

[ params_reduced, paramsP_reduced, reduced_rhs, ...
    reduced_stiffness_matrix] = galerkin_extension( params, paramsP, ...
    transformed_grid, rhs_online, stiffness_matrix_online, ...
    reduced_basis_matrix_B_pressure, reduced_basis_matrix_B_velocity);
velocity_reduced_dofs_online_parameter = ...
    reduced_basis_matrix_B_velocity * params_reduced.dofs;
pressure_reduced_dofs_online_parameter = ...
    reduced_basis_matrix_B_pressure * paramsP_reduced.dofs;

params_reduced_full = params;
params_reduced_full.dofs = velocity_reduced_dofs_online_parameter;
params_reduced_full = paramsP;
paramsP_reduced_full.dofs = pressure_reduced_dofs_online_parameter;

t_end_rb = toc();

% DG / full order solution for given parameter
tic();

% Creating stiffness matrix
[ stiffness_matrix_offline2, params] = assemble_stiffness_matrix_extension...
    ( params, paramsP, transformed_grid, mu, c11);

% Creating rhs vector
[ rhs_offline2, params] = assemble_rhs_extension...
    ( params, paramsP, transformed_grid, mu, c11);

close all

[ params, paramsP] = solve_plot_solution_schur...
    ( params, paramsP, transformed_grid, rhs_offline2, ...
    stiffness_matrix_offline2);
t_end_full = toc();
disp(['Time for offline computation : ' num2str(t_end)]);

velocity_dg_dofs_online_parameter = params.dofs;
pressure_dg_dofs_online_parameter = paramsP.dofs;


tic();
params_reduced_full = params;
params_reduced_full.dofs = velocity_reduced_dofs_online_parameter;
paramsP_reduced_full = paramsP;
paramsP_reduced_full.dofs = pressure_reduced_dofs_online_parameter;

for i=1:1:params_reduced_full.dimrange
    figure()
    axis equal
    [scalar_dofs, scalar_df_info] = ...
        ldg_scalar_component(params_reduced_full.dofs,i,params_reduced_full);
    sdf = ldgdiscfunc(scalar_dofs,scalar_df_info);
    disp(['Plotting ',num2str(i),' degree of freedom (Schur)'])
    %subplot(params.dimrange,1,i)
    %title(['Velocity degree of freedom number ',num2str(i)])
    if i==1
        title(['Velocity in x direction RB solution (Schur)'])
    else
        title(['Velocity in y direction RB solution (Schur)'])
    end
    axis equal
    axis tight
    ldg_plot(sdf,transformed_grid,params_reduced_full);
    plot(transformed_grid);
end

for i=1:1:paramsP_reduced_full.dimrange
    figure()
    [scalar_dofs, scalar_df_info] = ...
        ldg_scalar_component(paramsP_reduced_full.dofs,i,paramsP_reduced_full);
    sdf = ldgdiscfunc(scalar_dofs,scalar_df_info);
    disp(['Plotting ',num2str(i),' degree of freedom (for pressure)'])
    %subplot(paramsP.dimrange,1,i)
    title('Pressure RB solution (Schur)')
    axis equal
    axis tight
    ldg_plot(sdf,transformed_grid,paramsP_reduced_full);
    plot(transformed_grid);
end

t_end_rb_2 = toc();

disp(['Time for online computation : ' num2str(t_end_rb + t_end_rb_2)]);

error_rb_l2_velocity_online_parameter = error_rb_l2_norm_assembly...
    ( params, transformed_grid, velocity_reduced_dofs_online_parameter) / ...
    error_rb_l2_norm_assembly...
    ( params, transformed_grid, zeros(params.ndofs,1));
error_rb_l2_pressure_online_parameter = error_rb_l2_norm_assembly...
    ( paramsP, transformed_grid, pressure_reduced_dofs_online_parameter) / ...
    error_rb_l2_norm_assembly...
    ( paramsP, transformed_grid, zeros(paramsP.ndofs,1));

disp(['Speed up factor : ' num2str(t_end_full/(t_end_rb+t_end_rb_2))]);
disp(['RB l2 error for velocity : ' ...
    num2str(error_rb_l2_velocity_online_parameter)]);
disp(['RB l2 error for pressure : ' ...
    num2str(error_rb_l2_pressure_online_parameter)]);