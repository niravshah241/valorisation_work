clc
close all
clear all
%% grid making and saving
% pdetool
% draw geometry and save file as pdetool1_extension.m or similar
% pdetool1_extension.m or similar
% mesh --> export mesh
% save('my_grids_extension/mygridnirav1_extension','p','e','t')
%% start grid structure
% DO NOT DELETE : actual grid
params.mesh_number = 3;
params.gridtype = 'triagrid';
params.grid_show = false;
params.reference_parameter = [0.5 0.3];

tic();
%[ grid, params, el_subd] = grid_dd_func_extension(params);
%[ grid, params, el_subd ] = ...
%     grid_dd_func_extension_lid_driven_cavity( params);
%[ grid, params, el_subd] = grid_dd_func_extension_analytical(params);
%% Rectangular analytical grid
params.xnumintervals = 10;
params.ynumintervals = 10;
params.xrange = [0,1];
params.yrange = [0,1];
params.bnd_rect_corner1 = [-1,-1;1-eps,1/(params.ynumintervals*2)+1e-14]';
params.bnd_rect_corner2 = [2,2;1+eps,1-eps]';
params.bnd_rect_index = [-1,-2];
grid = construct_grid(params);

t_end = toc();
disp(['Time for grid generation : ' num2str(t_end)]);

nrep = [3 6 10 15];

params.pdeg = 2;
paramsP.pdeg = params.pdeg - 1;
params.dimrange = 2;
paramsP.dimrange = 1;
params.qdeg = 2;
paramsP.qdeg = params.qdeg;
params.mu = 1e-6;
paramsP.mu = params.mu;
params.nelements = grid.nelements;
paramsP.nelements = grid.nelements;

params.ndofs_per_element= nrep(params.pdeg) * params.dimrange;
params.ndofs = params.ndofs_per_element * grid.nelements;
params.dofs = zeros(params.ndofs,1);

paramsP.ndofs_per_element= nrep(paramsP.pdeg) * paramsP.dimrange;
paramsP.ndofs = paramsP.ndofs_per_element * grid.nelements;
paramsP.dofs = zeros(paramsP.ndofs,1);
params.show_sparsity = false;
paramsP.show_sparsity = params.show_sparsity;

params.rhs_func = @(glob,params,paramsP,grid)[ 2 * params.mu - 1 0]';

c11_range = [linspace(2e-4,2e-2,20) linspace(2e-2,1,10)];

for i = 1:1:length(c11_range)
    disp(['Entering in ' num2str(i) ' of ' num2str(length(c11_range))])
    c11 = c11_range(i); %1e-2;
    qdeg = params.qdeg;
    mu = params.mu;
    
    % check_assembly_n_u_tensor_jump_n_phi_tensor_jump_extension;
    % c11_calculation;
    % disp('TODO : check_assembly_n_u_tensor_jump_n_phi_tensor_jump_extension');
    % disp('TODO : c11_calculation');
    
    %% Creating stiffness matrix
    [ bilinear_side, diffusive_term, coercivity_term, flux_approximation] = ...
        assemble_bilinear_side_extension( params, paramsP, grid, mu, c11 );
    [ bilinear_side_pressure_terms, incompressibility_term, pressure_flux_term ] = ...
        assemble_bilinear_side_pressure_terms_extension...
        ( params, paramsP, grid);
    stiffness_matrix_offline = sparse(zeros(params.ndofs+paramsP.ndofs));
    stiffness_matrix_offline(1:params.ndofs,1:params.ndofs) = bilinear_side;
    stiffness_matrix_offline(1:params.ndofs,...
        params.ndofs+1:params.ndofs+paramsP.ndofs) = ...
        bilinear_side_pressure_terms';
    stiffness_matrix_offline(params.ndofs+1:params.ndofs+paramsP.ndofs,...
        1:params.ndofs) = bilinear_side_pressure_terms;
    
    
    %% Creating rhs vector
    rhs = zeros(params.ndofs + paramsP.ndofs,1);
    [ linear_side_offline, source_vector_offline ] = ...
        assemble_linear_side_extension( params, paramsP, grid, mu, c11 );
    linear_side_continuity = assemble_rhs_continuity( params, paramsP, grid);
    rhs(1:params.ndofs) = linear_side_offline;
    rhs(params.ndofs+1:params.ndofs+paramsP.ndofs) = linear_side_continuity;
    
    params.bilinear_side = bilinear_side;
    params.bilinear_side_pressure_terms = bilinear_side_pressure_terms';
    params.linear_side = rhs(1:params.ndofs);
    params.rhs_continuity = rhs(params.ndofs+1:params.ndofs+paramsP.ndofs);
    close all
    [ params, paramsP] = solve_plot_solution_schur...
        ( params, paramsP, grid, rhs, stiffness_matrix_offline);
    
    params.dof_analytical = @(glob) [glob(2) * (1 - glob(2)) 0];
    paramsP.dof_analytical = @(glob) (1 - glob(1));
    
    velocity_l2_error = error_l2_norm_assembly( params, grid);
    pressure_l2_error = error_l2_norm_assembly( paramsP, grid);
    matrix_condition_number = condest(stiffness_matrix_offline);
    
    velocity_error(i) = velocity_l2_error;
    pressure_error(i) = pressure_l2_error;
    condition_number(i) = matrix_condition_number;
end
close all
figure()
plot(c11_range,velocity_error,'*')
axis tight
xlabel('c11')
ylabel('velocity error')
title('c11 vs velocity error')

figure()
plot(c11_range,pressure_error,'*')
axis tight
xlabel('c11')
ylabel('pressure error')
title('c11 vs pressure error')

figure()
plot(c11_range,condition_number,'*')
axis tight
xlabel('c11')
ylabel('Condition number')
title('c11 vs condition number')

% % % Parametrization
% % N = 40;
% % x_para = 0.4 + (0.6-0.4).*rand(N,1);
% % y_para = 0.2 + (0.4-0.2).*rand(N,1);
% % snapshot_matrix_velocity = zeros(params.ndofs,N);
% % snapshot_matrix_pressure = zeros(paramsP.ndofs,N);
% % 
% % for temp = 1:1:N
% %     disp(['Entering parameter number ' num2str(temp) ' of ' num2str(N)])
% %     mu_x = x_para(temp);
% %     mu_y = y_para(temp);
% % 
% %     tic();
% %     [ transformed_grid, F_transformation_matrix, ...
% %         C_translation_vector] = transform_grid( params, grid, mu_x, mu_y);
% %     t_end = toc();
% %     close all
% %     disp(['Time for grid transformation : ' num2str(t_end)]);
% % 
% %     for i = 1:1:length(F_transformation_matrix)
% %         para_mapping{i} = inv(F_transformation_matrix{i});
% %     end
% % 
% %     assert(length(para_mapping) == length(el_subd),...
% %         'Number of subdomains and number of parametric mappings not same');
% % 
% %     %para_mapping is map from \Omega_{\mu} to \hat{\Omega}
% %     %i.e para_mapping is T and not F
% % 
% %     % Online stiffness matrix and rhs
% % 
% %     [ bilinear_side_online ] = assemble_bilinear_side_extension_online...
% %         ( params, paramsP, grid, diffusive_term, coercivity_term, ...
% %         flux_approximation, para_mapping, el_subd, transformed_grid, c11, mu);
% %     [ bilinear_side_pressure_terms_online ] = ...
% %         assemble_bilinear_side_pressure_terms_extension_online...
% %         ( params, paramsP, grid, incompressibility_term, ...
% %         pressure_flux_term, para_mapping, el_subd);
% %     stiffness_matrix_online = sparse(zeros(params.ndofs+paramsP.ndofs));
% %     stiffness_matrix_online(1:params.ndofs,1:params.ndofs) = ...
% %         bilinear_side_online;
% %     stiffness_matrix_online(1:params.ndofs,...
% %         params.ndofs+1:params.ndofs+paramsP.ndofs) = ...
% %         bilinear_side_pressure_terms_online;
% %     stiffness_matrix_online(params.ndofs+1:params.ndofs+paramsP.ndofs,...
% %         1:params.ndofs) = bilinear_side_pressure_terms_online';
% % 
% %     rhs_online = zeros(params.ndofs + paramsP.ndofs,1);
% %     linear_side_online = assemble_linear_side_extension_online...
% %         ( params, paramsP, grid, source_vector_offline, ...
% %         linear_side_offline, el_subd, para_mapping );
% %     rhs_online(1:params.ndofs,1) = linear_side_online;
% %     rhs_online(params.ndofs+1:params.ndofs+paramsP.ndofs,1) = ...
% %         linear_side_continuity;
% %     params.bilinear_side = bilinear_side_online;
% %     params.bilinear_side_pressure_terms = ...
% %         bilinear_side_pressure_terms_online;
% %     params.linear_side = rhs_online(1:params.ndofs);
% %     params.rhs_continuity = rhs_online(params.ndofs+1:...
% %         params.ndofs + paramsP.ndofs);
% % 
% %     [ params, paramsP] = solve_plot_solution_schur...
% %         ( params, paramsP, transformed_grid, rhs_online, ...
% %         stiffness_matrix_online);
% % 
% %     snapshot_matrix_velocity(:,temp) = params.dofs;
% %     snapshot_matrix_pressure(:,temp) = paramsP.dofs;
% % 
% % end
% % 
% % % POD-Galerkin
% % 
% % min_eigen_value = 0;
% % max_reduced_basis = 30;
% % inner_product_matrix = ldg_mass_matrix(params,grid,params);
% % 
% % [ reduced_basis_matrix_B_velocity, eigen_values_velocity] = ...
% %     pod_extension( snapshot_matrix_velocity, params, min_eigen_value, ...
% %     max_reduced_basis, inner_product_matrix );
% % 
% % min_eigen_value = 0;
% % max_reduced_basis = 3;
% % inner_product_matrix = ldg_mass_matrix(paramsP,grid,paramsP);
% % 
% % [ reduced_basis_matrix_B_pressure, eigen_values_pressure] = ...
% %     pod_extension( snapshot_matrix_pressure, paramsP, min_eigen_value, ...
% %     max_reduced_basis, inner_product_matrix );
% % 
% % % Galerkin projection and rb error
% % 
% % N = 10;
% % x_para = 0.4 + (0.6-0.4).*rand(N,1);
% % y_para = 0.2 + (0.4-0.2).*rand(N,1);
% % error_rb_l2_velocity = zeros(N,1);
% % error_rb_l2_pressure = zeros(N,1);
% % 
% % for temp = 1:1:N
% %     disp(['Entering parameter number ' num2str(temp) ' of ' num2str(N)])
% %     mu_x = x_para(temp);
% %     mu_y = y_para(temp);
% % 
% %     tic();
% %     [ transformed_grid, F_transformation_matrix, ...
% %         C_translation_vector] = transform_grid( params, grid, mu_x, mu_y);
% %     t_end = toc();
% %     close all
% %     disp(['Time for grid transformation : ' num2str(t_end)]);
% % 
% %     for i = 1:1:length(F_transformation_matrix)
% %         para_mapping{i} = inv(F_transformation_matrix{i});
% %     end
% % 
% %     assert(length(para_mapping) == length(el_subd),...
% %         'Number of subdomains and number of parametric mappings not same');
% % 
% %     para_mapping is map from \Omega_{\mu} to \hat{\Omega}
% %     i.e para_mapping is T and not F
% % 
% %     % Online stiffness matrix and rhs
% % 
% %     [ bilinear_side_online ] = assemble_bilinear_side_extension_online...
% %         ( params, paramsP, grid, diffusive_term, coercivity_term, ...
% %         flux_approximation, para_mapping, el_subd, transformed_grid, ...
% %         c11, mu);
% %     [ bilinear_side_pressure_terms_online ] = ...
% %         assemble_bilinear_side_pressure_terms_extension_online...
% %         ( params, paramsP, grid, incompressibility_term, ...
% %         pressure_flux_term, para_mapping, el_subd);
% %     stiffness_matrix_online = sparse(zeros(params.ndofs+paramsP.ndofs));
% %     stiffness_matrix_online(1:params.ndofs,1:params.ndofs) = ...
% %         bilinear_side_online;
% %     stiffness_matrix_online(1:params.ndofs,...
% %         params.ndofs+1:params.ndofs+paramsP.ndofs) = ...
% %         bilinear_side_pressure_terms_online;
% %     stiffness_matrix_online(params.ndofs+1:params.ndofs+paramsP.ndofs,...
% %         1:params.ndofs) = bilinear_side_pressure_terms_online';
% % 
% %     rhs_online = zeros(params.ndofs + paramsP.ndofs,1);
% %     linear_side_online = assemble_linear_side_extension_online...
% %         ( params, paramsP, grid, source_vector_offline, ...
% %         linear_side_offline, el_subd, para_mapping );
% %     rhs_online(1:params.ndofs,1) = linear_side_online;
% %     rhs_online(params.ndofs+1:params.ndofs+paramsP.ndofs,1) = ...
% %         linear_side_continuity;
% %     params.bilinear_side = bilinear_side_online;
% %     params.bilinear_side_pressure_terms = ...
% %         bilinear_side_pressure_terms_online;
% %     params.linear_side = rhs_online(1:params.ndofs);
% %     params.rhs_continuity = rhs_online(params.ndofs+1:...
% %         params.ndofs + paramsP.ndofs);
% % 
% %     [ params, paramsP] = solve_plot_solution_schur...
% %         ( params, paramsP, transformed_grid, rhs_online, ...
% %         stiffness_matrix_online);
% % 
% %     [ params_reduced, paramsP_reduced, reduced_rhs, ...
% %         reduced_stiffness_matrix] = galerkin_extension( params, paramsP, ...
% %         transformed_grid, rhs_online, stiffness_matrix_online, ...
% %         reduced_basis_matrix_B_pressure, reduced_basis_matrix_B_velocity);
% %     velocity_dofs = reduced_basis_matrix_B_velocity * params_reduced.dofs;
% %     pressure_dofs = reduced_basis_matrix_B_pressure * paramsP_reduced.dofs;
% %     error_rb_l2_velocity(temp) = error_rb_l2_norm_assembly...
% %         ( params, transformed_grid, velocity_dofs);
% %     error_rb_l2_pressure(temp) = error_rb_l2_norm_assembly...
% %         ( paramsP, transformed_grid, pressure_dofs);
% % end
% % 
% % error_rb_l2_velocity_mean = mean(error_rb_l2_velocity);
% % error_rb_l2_pressure_mean = mean(error_rb_l2_pressure);
% % 
% % % Stiffness matrix and rhs check
% % [ params, paramsP, rhs_real, stifness_matrix_real] = ...
% %     assemble_stifness_matrix...
% %     ( params, paramsP, grid, params.qdeg, mu, c11 );
% % error_rhs_assembly = max(abs(rhs - rhs_real));
% % error_stiffness_matrix_assembly = max(max(abs(stiffness_matrix_offline ...
% %     - stifness_matrix_real)));
% % 
% % [ params, paramsP, rhs_real, stifness_matrix_real] = ...
% %     assemble_stifness_matrix...
% %     ( params, paramsP, transformed_grid, params.qdeg, mu, c11 );
% % error_rhs_affine = max(abs(rhs_online - rhs_real));
% % error_stiffness_matrix_affine = max(max(abs(stiffness_matrix_online ...
% %     - stifness_matrix_real)));
% % 
% % [ params, paramsP, rhs_real, stifness_matrix_real] = ...
% %     assemble_stifness_matrix...
% %     ( params, paramsP, grid, params.qdeg, mu, c11 );
% % error_rhs_assembly = max(abs(rhs - rhs_real));
% % error_stiffness_matrix_assembly = max(max(abs(stiffness_matrix_offline ...
% %     - stifness_matrix_real)));
% % 
% % [ params, paramsP, rhs_real, stifness_matrix_real] = ...
% %     assemble_stifness_matrix...
% %     ( params, paramsP, transformed_grid, params.qdeg, mu, c11 );
% % error_rhs_affine = max(abs(rhs_online - rhs_real));
% % error_stiffness_matrix_affine = max(max(abs(stiffness_matrix_online ...
% %     - stifness_matrix_real)));
% % 
% % [ params, paramsP] = solve_plot_solution_schur...
% %     ( params, paramsP, grid, rhs, stiffness_matrix_offline);