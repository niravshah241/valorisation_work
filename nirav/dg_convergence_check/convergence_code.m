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
num_intervals = [4 6 8 10 12 14];
relative_error_l2_velocity = zeros(length(num_intervals),1);
relative_error_l2_pressure = zeros(length(num_intervals),1);
for i = 1:1:length(num_intervals)
    disp(['Step size : ' num2str(num_intervals(i))]);
    tic();
    params.gridtype = 'triagrid';
    params.grid_show = false;
    
    tic();
    % [ grid, params, el_subd] = grid_dd_func_extension(params);
    % [ grid, params, el_subd ] = ...
    %     grid_dd_func_extension_lid_driven_cavity( params);
    %% Rectangular analytical grid
    params.xnumintervals = num_intervals(i);
    params.ynumintervals = num_intervals(i);
    params.xrange = [0,1];
    params.yrange = [0,1];
    params.bnd_rect_corner1 = [-1,-1;1-eps,1/(params.ynumintervals*2)+1e-14]';
    params.bnd_rect_corner2 = [2,2;1+eps,1-eps]';
    params.bnd_rect_index = [-1,-2];
    grid = construct_grid(params);
    % params.grid_initfile = ['convergence_mesh_',...
    %     num2str(params.mesh_number),'.mat'];
    % load(params.grid_initfile)
    % params.bnd_rect_corner1 = [-1,-1;1-eps,eps]';%1/(params.ynumintervals*2)+1e-14]';
    % params.bnd_rect_corner2 = [2,2;1+eps,1-eps]';
    % params.bnd_rect_index = [-1,-2];
    % grid = construct_grid(params);
    
    t_end = toc();
    disp(['Time for grid generation : ' num2str(t_end)]);
    
    nrep = [3 6 10 15];
    
    params.pdeg = 2;
    paramsP.pdeg = params.pdeg - 1;
    params.dimrange = 2;
    paramsP.dimrange = 1;
    params.qdeg = 2;
    paramsP.qdeg = params.qdeg;
    params.mu = 1e-3;
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
    
    params.rhs_func = @(glob,params,paramsP,grid)[2*params.mu-1 0]';
    
    tic();
    c11 = 1e4;
    qdeg = params.qdeg;
    mu = params.mu;
    
    % Creating stiffness matrix
    tic();
    [ stiffness_matrix_offline, params, bilinear_side, diffusive_term, ...
        coercivity_term, flux_approximation, bilinear_side_pressure_terms, ...
        incompressibility_term, pressure_flux_term ] = ...
        assemble_stiffness_matrix_extension...
        ( params, paramsP, grid, mu, c11);
    
    % Creating rhs vector
    [ rhs_offline, params, linear_side_offline, source_vector_offline, ...
        linear_side_continuity ] = assemble_rhs_extension...
        ( params, paramsP, grid, mu, c11);
    assembly_time = toc();
    disp(['Assembly time of lhs matrix and rhs vector : ' num2str(assembly_time)]);
    
    close all
    tic();
    [ params, paramsP] = solve_plot_solution_schur...
        ( params, paramsP, grid, rhs_offline, stiffness_matrix_offline);
    solution_time = toc();
    disp(['Time for solving system : ' num2str(solution_time)]);
    
    params.dof_analytical = @(glob) (glob(2) * (1 - glob(2)));
    paramsP.dof_analytical = @(glob) (1 - glob(1));
    velocity_dofs = params.dofs;
    pressure_dofs = paramsP.dofs;
    
    % Computing velocity error
    absolute_error_l2_velocity = error_l2_norm_assembly(params,grid);
    params.dofs = zeros(params.ndofs,1);
    relative_error_l2_velocity(i) = absolute_error_l2_velocity / ...
        error_l2_norm_assembly(params,grid);
    params.dofs = velocity_dofs;
    
    % Computing pressure error
    absolute_error_l2_pressure = error_l2_norm_assembly(paramsP,grid);
    paramsP.dofs = zeros(paramsP.ndofs,1);
    relative_error_l2_pressure(i) = absolute_error_l2_pressure / ...
        error_l2_norm_assembly(paramsP,grid);
    paramsP.dofs = pressure_dofs;
    
    total_run_time = toc();
    disp(['Total run time : ' num2str(total_run_time)]);
end

a = figure();
semilogy(1./num_intervals,relative_error_l2_velocity,'*-');
axis tight
xlabel('Step size');
ylabel('Relative error L^2 velocity (log scale)');
str = 'Relative error L^2 velocity vs step size (semilog scale)),';
title([ str, newline, '\nu = ',num2str(params.mu),' c_{11} = ',num2str(c11)]);
saveas(a,'nirav/dg_convergence_check/step_size_vs_velocity_l2_error_semilog.fig');
saveas(a,'nirav/dg_convergence_check/step_size_vs_velocity_l2_error_semilog.jpg');

a = figure();
semilogy(1./num_intervals,relative_error_l2_pressure,'*-');
axis tight
xlabel('Step size');
ylabel('Relative error L^2 pressure (log scale)');
str = 'Relative error L^2 pressure vs step size (semilog scale)),';
title([ str, newline, '\nu = ',num2str(params.mu),' c_{11} = ',num2str(c11)]);
saveas(a,'nirav/dg_convergence_check/step_size_vs_pressure_l2_error_semilog.fig');
saveas(a,'nirav/dg_convergence_check/step_size_vs_pressure_l2_error_semilog.jpg');

a = figure();
loglog(1./num_intervals,relative_error_l2_velocity,'*-');
axis tight
xlabel('Step size (log scale)');
ylabel('Relative error L^2 velocity (log scale)');
str = 'Relative error L^2 velocity vs step size (loglog scale)),';
title([ str, newline, '\nu = ',num2str(params.mu),' c_{11} = ',num2str(c11)]);
saveas(a,'nirav/dg_convergence_check/step_size_vs_velocity_l2_error_loglog.fig');
saveas(a,'nirav/dg_convergence_check/step_size_vs_velocity_l2_error_loglog.jpg');

a = figure();
loglog(1./num_intervals,relative_error_l2_pressure,'*-');
axis tight
xlabel('Step size (log scale)');
ylabel('Relative error L^2 pressure (log scale)');
str = 'Relative error L^2 pressure vs step size (loglog scale)),';
title([ str, newline, '\nu = ',num2str(params.mu),' c_{11} = ',num2str(c11)]);
saveas(a,'nirav/dg_convergence_check/step_size_vs_pressure_l2_error_loglog.fig');
saveas(a,'nirav/dg_convergence_check/step_size_vs_pressure_l2_error_loglog.jpg');

save('nirav/dg_convergence_check/relative_error_l2_velocity',...
    'relative_error_l2_velocity');
save('nirav/dg_convergence_check/relative_error_l2_pressure',...
    'relative_error_l2_pressure');
step_size = 1./num_intervals;
save('nirav/dg_convergence_check/step_size','step_size');