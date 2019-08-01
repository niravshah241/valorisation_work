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
[ grid, params, el_subd] = grid_dd_func_extension(params);
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

params.rhs_func = @(glob,params,paramsP,grid)[12 10]';

c11 = 1e-2;
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

%% Parametrization
N = 10;
x_para = 0.48 + (0.52-0.48).*rand(N,1);
y_para = 0.28 + (0.32-0.28).*rand(N,1);
snapshot_matrix_velocity = zeros(params.ndofs,N);
snapshot_matrix_pressure = zeros(paramsP.ndofs,N);

for temp = 1:1:N
    disp(['Entering parameter number ' num2str(temp) ' of ' num2str(N)])
    mu_x = x_para(temp);
    mu_y = y_para(temp);
    
    tic();
    [ transformed_grid, F_transformation_matrix, ...
        C_translation_vector] = transform_grid( params, grid, mu_x, mu_y);
    t_end = toc();
    close all
    disp(['Time for grid transformation : ' num2str(t_end)]);
    
    for i = 1:1:length(F_transformation_matrix)
        para_mapping{i} = inv(F_transformation_matrix{i});
    end
    
    assert(length(para_mapping) == length(el_subd),...
        'Number of subdomains and number of parametric mappings not same');
    
    %para_mapping is map from \Omega_{\mu} to \hat{\Omega}
    %i.e para_mapping is T and not F
    
    %% Online stiffness matrix and rhs
    
    [ bilinear_side_online ] = assemble_bilinear_side_extension_online...
        ( params, paramsP, grid, diffusive_term, coercivity_term, ...
        flux_approximation, para_mapping, el_subd, transformed_grid, c11, mu);
    [ bilinear_side_pressure_terms_online ] = ...
        assemble_bilinear_side_pressure_terms_extension_online...
        ( params, paramsP, grid, incompressibility_term, ...
        pressure_flux_term, para_mapping, el_subd);
    stiffness_matrix_online = sparse(zeros(params.ndofs+paramsP.ndofs));
    stiffness_matrix_online(1:params.ndofs,1:params.ndofs) = ...
        bilinear_side_online;
    stiffness_matrix_online(1:params.ndofs,...
        params.ndofs+1:params.ndofs+paramsP.ndofs) = ...
        bilinear_side_pressure_terms_online;
    stiffness_matrix_online(params.ndofs+1:params.ndofs+paramsP.ndofs,...
        1:params.ndofs) = bilinear_side_pressure_terms_online';
    
    rhs_online = zeros(params.ndofs + paramsP.ndofs,1);
    linear_side_online = assemble_linear_side_extension_online...
        ( params, paramsP, grid, source_vector_offline, ...
        linear_side_offline, el_subd, para_mapping );
    rhs_online(1:params.ndofs,1) = linear_side_online;
    rhs_online(params.ndofs+1:params.ndofs+paramsP.ndofs,1) = ...
        linear_side_continuity;
    params.bilinear_side = bilinear_side_online;
    params.bilinear_side_pressure_terms = ...
        bilinear_side_pressure_terms_online;
    params.linear_side = rhs_online(1:params.ndofs);
    params.rhs_continuity = rhs_online(params.ndofs+1:...
        params.ndofs + paramsP.ndofs);
    
    [ params, paramsP] = solve_plot_solution_schur...
        ( params, paramsP, transformed_grid, rhs_online, ...
        stiffness_matrix_online);
    
    snapshot_matrix_velocity(:,temp) = params.ndofs;
    snapshot_matrix_pressure(:,temp) = paramsP.ndofs;
    
end

%% POD-Galerkin

min_eigen_value = 0;
max_reduced_basis = 3;
inner_product_matrix = ldg_mass_matrix(params,grid,params);

[ reduced_basis_matrix_B_velocity, eigen_values_velocity] = ...
    pod_extension( snapshot_matrix_velocity, params, min_eigen_value, ...
    max_reduced_basis, inner_product_matrix );

min_eigen_value = 0;
max_reduced_basis = 3;
inner_product_matrix = ldg_mass_matrix(paramsP,grid,paramsP);

[ reduced_basis_matrix_B_pressure, eigen_values_pressure] = ...
    pod_extension( snapshot_matrix_pressure, paramsP, min_eigen_value, ...
    max_reduced_basis, inner_product_matrix );

%% Galerkin projection


%% Stiffness matrix and rhs check
% [ params, paramsP, rhs_real, stifness_matrix_real] = ...
%     assemble_stifness_matrix...
%     ( params, paramsP, grid, params.qdeg, mu, c11 );
% error_rhs_assembly = max(abs(rhs - rhs_real));
% error_stiffness_matrix_assembly = max(max(abs(stiffness_matrix_offline ...
%     - stifness_matrix_real)));
%
% [ params, paramsP, rhs_real, stifness_matrix_real] = ...
%     assemble_stifness_matrix...
%     ( params, paramsP, transformed_grid, params.qdeg, mu, c11 );
% error_rhs_affine = max(abs(rhs_online - rhs_real));
% error_stiffness_matrix_affine = max(max(abs(stiffness_matrix_online ...
%     - stifness_matrix_real)));
%
% [ params, paramsP, rhs_real, stifness_matrix_real] = ...
%     assemble_stifness_matrix...
%     ( params, paramsP, grid, params.qdeg, mu, c11 );
% error_rhs_assembly = max(abs(rhs - rhs_real));
% error_stiffness_matrix_assembly = max(max(abs(stiffness_matrix_offline ...
%     - stifness_matrix_real)));
%
% [ params, paramsP, rhs_real, stifness_matrix_real] = ...
%     assemble_stifness_matrix...
%     ( params, paramsP, transformed_grid, params.qdeg, mu, c11 );
% error_rhs_affine = max(abs(rhs_online - rhs_real));
% error_stiffness_matrix_affine = max(max(abs(stiffness_matrix_online ...
%     - stifness_matrix_real)));
%
% [ params, paramsP] = solve_plot_solution_schur...
%     ( params, paramsP, grid, rhs, stiffness_matrix_offline);