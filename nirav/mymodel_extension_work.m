clc
close all
clear all
clock_start_time = clock;
total_start_time = tic();
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
% [ grid, params, el_subd ] = ...
%     grid_dd_func_extension_lid_driven_cavity( params);

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
params.rhs_func = @(glob,params,paramsP,grid)[0 0]';

tic();
c11 = 1e-2;
qdeg = params.qdeg;
mu = params.mu;

tic();
% Creating stiffness matrix
[ stiffness_matrix_offline, params, bilinear_side, diffusive_term, ...
    coercivity_term, flux_approximation, bilinear_side_pressure_terms, ...
    incompressibility_term, pressure_flux_term ] = ...
    assemble_stiffness_matrix_extension...
    ( params, paramsP, grid, mu, c11);

% Creating rhs vector

[ rhs_offline, params, linear_side_offline, source_vector_offline, ...
    linear_side_continuity ] = assemble_rhs_extension...
    ( params, paramsP, grid, mu, c11);
close all
t_DG_assembly = toc();
disp(['Time for assembling DG stiffness matrix and rhs : ' num2str(t_DG_assembly)]);

tic();
[ params, paramsP] = solve_plot_solution_schur...
    ( params, paramsP, grid, rhs_offline, stiffness_matrix_offline);
t_DG_solve = toc();
disp(['Time for solving DG system : ' num2str(t_DG_solve)]);
disp(['Time for assembling+solving DG system : ' ...
    num2str(t_DG_assembly+t_DG_solve)]);

%%%% assembling offline stifness matrix and rhs
%%%% This is true offline time as extension files are used for parametrization
%%%% purpose whereas real assemblies can be assembled faster,
%%%% if parametritaion is not needed.
%%%% The difference for assembly could be upto 3 times.
% tic();
% [ params, paramsP, rhs, stifness_matrix] = ...
%     assemble_stifness_matrix( params, paramsP, grid, params.qdeg, mu, c11 );
% t_DG_assembly_short = toc();
% disp(['Time for assembling DG stiffness matrix and rhs : ' num2str(t_DG_assembly_short)]);

% % Parametrization
N = 20;
x_para = 0.4 + (0.6-0.4).*rand(N,1);
y_para = 0.2 + (0.4-0.2).*rand(N,1);
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
    
    % Online stiffness matrix and rhs
    [ stiffness_matrix_online, params ] = ...
        assemble_stiffness_matrix_online_extension( params, paramsP, ...
        grid, diffusive_term, coercivity_term, ...
        flux_approximation, para_mapping, el_subd, transformed_grid, c11, ...
        mu, incompressibility_term, pressure_flux_term);
    
    [ rhs_online, params ] = assemble_rhs_online_extension...
        ( params, paramsP, grid, source_vector_offline, ...
        linear_side_offline, el_subd, para_mapping, linear_side_continuity);
    
    [ params, paramsP] = solve_plot_solution_schur...
        ( params, paramsP, transformed_grid, rhs_online, ...
        stiffness_matrix_online);
    
    snapshot_matrix_velocity(:,temp) = params.dofs;
    snapshot_matrix_pressure(:,temp) = paramsP.dofs;
end

% POD-Galerkin
k = 1:1:min(size(snapshot_matrix_velocity,2),20);
error_rb_velocity_mean = zeros(length(k),1);
error_rb_pressure_mean = zeros(length(k),1);
error_rb_velocity_max = zeros(length(k),1);
error_rb_pressure_max = zeros(length(k),1);
error_rb_energy = zeros(length(k),1);
% error_estimate_velocity = zeros(length(k),1);
% error_estimate_pressure = zeros(length(k),1);
speedup = zeros(length(k),1);
online_simulation_time = zeros(length(k),1);

for temp2 = 1:1:length(k)
    disp(['Entering POD-Galerkin number ' num2str(temp2) ' of ' num2str(length(k))])
    min_eigen_value = [0 0];
    max_reduced_basis = [k(temp2) k(temp2)];
    
    inner_product_matrix_velocity = ldg_mass_matrix(params,grid,params) ...
        + ldg_h1_seminorm_mass_matrix_assembly(params,grid);
    
    [ reduced_basis_matrix_B_velocity, eigen_values_velocity, ...
        ~] = pod_extension_velocity...
        ( snapshot_matrix_velocity, params, min_eigen_value, ...
        max_reduced_basis, inner_product_matrix_velocity );
    
    %     disp(['Max error estimate for velocity : ' num2str(error_estimate_velocity{temp2})]);
    %     disp(['Before comparing error, please verify that ' newline...
    %         'error estimator and actual error are measured in same norm'])
    
    min_eigen_value = 0;
    max_reduced_basis = k(temp2);
    
    inner_product_matrix_pressure = ldg_mass_matrix(paramsP,grid,paramsP);
    
    [ reduced_basis_matrix_B_pressure, eigen_values_pressure, ...
        ~] = ...
        pod_extension( snapshot_matrix_pressure, paramsP, min_eigen_value, ...
        max_reduced_basis, inner_product_matrix_pressure);
    
    %     disp(['Max error estimate for pressure : ' num2str(error_estimate_pressure(temp2))]);
    %     disp(['Before comparing error, please verify that ' newline...
    %         'error estimator and actual error are measured in same norm'])
    
    % Galerkin projection and rb error
    
    N = 10;
    x_para = 0.4 + (0.6-0.4).*rand(N,1);
    y_para = 0.2 + (0.4-0.2).*rand(N,1);
    error_rb_velocity = zeros(N,1);
    error_rb_pressure = zeros(N,1);
    
    for temp = 1:1:N
        disp(['In POD-Galerkin number ' num2str(temp2) ' of ' num2str(length(k))])
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
        
        % %     para_mapping is map from \Omega_{\mu} to \hat{\Omega}
        % %     i.e para_mapping is T and not F
        
        % Online stiffness matrix and rhs
        
        [ stiffness_matrix_online, params ] = ...
            assemble_stiffness_matrix_online_extension( params, paramsP, ...
            grid, diffusive_term, coercivity_term, ...
            flux_approximation, para_mapping, el_subd, transformed_grid, c11, ...
            mu, incompressibility_term, pressure_flux_term);
        
        [ rhs_online, params ] = assemble_rhs_online_extension...
            ( params, paramsP, grid, source_vector_offline, ...
            linear_side_offline, el_subd, para_mapping, linear_side_continuity);
        
        [ params, paramsP] = solve_plot_solution_schur...
            ( params, paramsP, transformed_grid, rhs_online, ...
            stiffness_matrix_online);
        
        [ params_reduced, paramsP_reduced, reduced_rhs, ...
            reduced_stiffness_matrix] = galerkin_extension( params, paramsP, ...
            transformed_grid, rhs_online, stiffness_matrix_online, ...
            reduced_basis_matrix_B_pressure, reduced_basis_matrix_B_velocity);
        velocity_dofs = reduced_basis_matrix_B_velocity * params_reduced.dofs;
        pressure_dofs = reduced_basis_matrix_B_pressure * paramsP_reduced.dofs;
        error_vector_velocity = params.dofs - velocity_dofs;
        error_vector_pressure = paramsP.dofs - pressure_dofs;
        error_rb_velocity(temp) = ...
            ((params.dofs - velocity_dofs)' * inner_product_matrix_velocity * ...
            (params.dofs - velocity_dofs)) / ...
            (params.dofs' * inner_product_matrix_velocity * params.dofs);
        error_rb_pressure(temp) = ...
            ((paramsP.dofs - pressure_dofs)' * inner_product_matrix_pressure * ...
            (paramsP.dofs - pressure_dofs)) / ...
            (paramsP.dofs' * inner_product_matrix_pressure * paramsP.dofs);
        error_rb_energy(temp) = ...
            (([error_vector_velocity;error_vector_pressure])' * ...
            stiffness_matrix_online * ...
            ([error_vector_velocity;error_vector_pressure])) / ...
            (([params.dofs;paramsP.dofs])' * ...
            stiffness_matrix_online * ([params.dofs;paramsP.dofs]));
    end
    
    error_rb_velocity_mean(temp2) = mean(error_rb_velocity);
    error_rb_pressure_mean(temp2) = mean(error_rb_pressure);
    error_rb_velocity_max(temp2) = max(error_rb_velocity);
    error_rb_pressure_max(temp2) = max(error_rb_pressure);
    error_rb_energy_mean(temp2) = mean(error_rb_energy);
    error_rb_energy_max(temp2) = max(error_rb_energy);
    
    a = figure();
    semilogy(real(eigen_values_velocity{1}),'*');
    axis tight
    xlabel('number of eigen value');
    ylabel('eigenvalues x-velocity (log scale)');
    title('x-velocity eigenvalues (semilog scale))');
    saveas(a,'nirav/pod_galerkin/x_velocity_eigen_value_semilog.fig');
    saveas(a,'nirav/pod_galerkin/x_velocity_eigen_value_semilog.jpg');
    
    a = figure();
    semilogy(real(eigen_values_velocity{2}),'*');
    axis tight
    xlabel('number of eigen value');
    ylabel('eigenvalues y-velocity (log scale)');
    title('y-velocity eigenvalues (semilog scale))');
    saveas(a,'nirav/pod_galerkin/y_velocity_eigen_value_semilog.fig');
    saveas(a,'nirav/pod_galerkin/y_velocity_eigen_value_semilog.jpg');
    
    
    a = figure();
    semilogy(real(eigen_values_pressure),'*');
    axis tight
    xlabel('number of eigen value');
    ylabel('eigenvalues pressure (log scale)');
    title('pressure eigenvalues (semilog scale))');
    saveas(a,'nirav/pod_galerkin/pressure_eigen_value_semilog.fig');
    saveas(a,'nirav/pod_galerkin/pressure_eigen_value_semilog.jpg');
    
    a = figure();
    scatter(x_para,y_para,[],error_rb_velocity,'filled');
    axis tight
    xlabel('x-coordinate of tip');
    ylabel('y-coordinate of tip');
    title('Velocity RB relative error');
    c = colorbar('eastoutside');
    c.Label.String='RB Velocity error';
    saveas(a,'nirav/pod_galerkin/rb_error_velocity.fig');
    saveas(a,'nirav/pod_galerkin/rb_error_velocity.jpg');
    
    a = figure();
    scatter(x_para,y_para,[],error_rb_pressure,'filled');
    axis tight
    xlabel('x-coordinate of tip');
    ylabel('y-coordinate of tip');
    title('Pressure RB relative error');
    c = colorbar('eastoutside');
    c.Label.String='RB Pressure error';
    saveas(a,'nirav/pod_galerkin/rb_error_pressure.fig');
    saveas(a,'nirav/pod_galerkin/rb_error_pressure.jpg');
    
    run_time = toc();
    disp(['Run time : ' num2str(run_time)]);
    
    %% %%%online phase
    mu_x = 0.43; % online parameter 1
    mu_y = 0.36; % online parameter 2
    plot_online_solution = 0;
    online_phase;
    speedup(temp2) = t_end_full/(t_end_rb+t_end_rb_2);
    online_simulation_time(temp2) = t_end_rb+t_end_rb_2;
end

a = figure();
semilogy(k,speedup,'*-');
axis tight
xlabel('Size of reduced basis');
ylabel('Speed up (log scale)');
title('Size of reduced basis vs Speed up (semilog scale)');
saveas(a,'nirav/pod_galerkin/size_vs_reduced_basis_semilog.fig');
saveas(a,'nirav/pod_galerkin/size_vs_reduced_basis_semilog.jpg');

a = figure();
semilogy(k,online_simulation_time,'*-');
axis tight
xlabel('Size of reduced basis');
ylabel('Online simulation time (log scale)');
title('Size of reduced basis vs Online simulation time (semilog scale)');
saveas(a,'nirav/pod_galerkin/size_vs_online_simulation_time_semilog.fig');
saveas(a,'nirav/pod_galerkin/size_vs_online_simulation_time_semilog.jpg');

a = figure();
semilogy(k,error_rb_velocity_max,'*-');
axis tight
xlabel('Size of reduced basis');
ylabel('Maximum reduced basis velocity error (log scale)');
title('Size of reduced basis vs maximum velocity error (semilog scale)');
saveas(a,'nirav/pod_galerkin/size_vs_maximum_reduced_basis_velocity_error_semilog.fig');
saveas(a,'nirav/pod_galerkin/size_vs_maximum_reduced_basis_velocity_error_semilog.jpg');

a = figure();
semilogy(k,error_rb_velocity_mean,'*-');
axis tight
xlabel('Size of reduced basis');
ylabel('Average reduced basis velocity error (log scale)');
title('Size of reduced basis vs average velocity error (semilog scale)');
saveas(a,'nirav/pod_galerkin/size_vs_average_reduced_basis_velocity_error_semilog.fig');
saveas(a,'nirav/pod_galerkin/size_vs_average_reduced_basis_velocity_error_semilog.jpg');

a = figure();
semilogy(k,error_rb_pressure_mean,'*-');
axis tight
xlabel('Size of reduced basis');
ylabel('Average reduced basis pressure error (log scale)');
title('Size of reduced basis vs average pressure error (semilog scale)');
saveas(a,'nirav/pod_galerkin/size_vs_average_reduced_basis_pressure_error_semilog.fig');
saveas(a,'nirav/pod_galerkin/size_vs_average_reduced_basis_pressure_error_semilog.jpg');

a = figure();
semilogy(k,error_rb_pressure_max,'*-');
axis tight
xlabel('Size of reduced basis');
ylabel('Maximum reduced basis pressure error (log scale)');
title('Size of reduced basis vs Maximum pressure error (semilog scale)');
saveas(a,'nirav/pod_galerkin/size_vs_maximum_reduced_basis_pressure_error_semilog.fig');
saveas(a,'nirav/pod_galerkin/size_vs_maximum_reduced_basis_pressure_error_semilog.jpg');

% a = figure();
% semilogy(k,error_rb_energy_mean,'*-');
% axis tight
% xlabel('Size of reduced basis');
% ylabel('Average reduced basis energy error (log scale)');
% title('Size of reduced basis vs average energy error (semilog scale)');
% saveas(a,'nirav/pod_galerkin/size_vs_average_reduced_basis_energy_error_semilog.fig');
% saveas(a,'nirav/pod_galerkin/size_vs_average_reduced_basis_energy_error_semilog.jpg');
%
% a = figure();
% semilogy(k,error_rb_energy_max,'*-');
% axis tight
% xlabel('Size of reduced basis');
% ylabel('Maximum reduced basis energy error (log scale)');
% title('Size of reduced basis vs Maximum energy error (semilog scale)');
% saveas(a,'nirav/pod_galerkin/size_vs_maximum_reduced_basis_energy_error_semilog.fig');
% saveas(a,'nirav/pod_galerkin/size_vs_maximum_reduced_basis_energy_error_semilog.jpg');

a = figure();
M = 50;
x_para_i = linspace(min(x_para),max(x_para),M);
y_para_i = linspace(min(y_para),max(y_para),M);
[Xi,Yi] = meshgrid(x_para_i,y_para_i);
Ci = griddata(x_para,y_para,error_rb_velocity,Xi,Yi);
mesh(Xi,Yi,Ci)
hold on
plot3(x_para,y_para,error_rb_velocity,'o')
xlabel('x-coordinate of tip')
ylabel('y-coordinate of tip')
zlabel('Relative error rb Velocity')
title('RB relative error velocity over parameter space')
saveas(a,'nirav/pod_galerkin/rb_error_velocity.fig');
saveas(a,'nirav/pod_galerkin/rb_error_velocity.jpg');

a = figure();
M = 50;
x_para_i = linspace(min(x_para),max(x_para),M);
y_para_i = linspace(min(y_para),max(y_para),M);
[Xi,Yi] = meshgrid(x_para_i,y_para_i);
Ci = griddata(x_para,y_para,error_rb_pressure,Xi,Yi);
mesh(Xi,Yi,Ci)
hold on
plot3(x_para,y_para,error_rb_pressure,'o')
xlabel('x-coordinate of tip')
ylabel('y-coordinate of tip')
zlabel('Relative error rb Pressure')
title('RB relative error pressure over parameter space')
saveas(a,'nirav/pod_galerkin/rb_error_pressure.fig');
saveas(a,'nirav/pod_galerkin/rb_error_pressure.jpg');

total_end_time = toc(total_start_time);
disp(['Total time for script : ' num2str(total_end_time)]);
etime(clock,clock_start_time)