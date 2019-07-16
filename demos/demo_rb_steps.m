% small collection of RBmatlab's rb-simulation abilities
% please inspect the source-file during execution, as the main
% purpose is to demonstrate the possibilities of RBmatlab's
% rb-philosophy

% Bernard Haasdonk 11.4.2007

help demo_rb_steps;

disp('Small demonstration of RBmatlab reduced basis commands');
disp(['Please type dbcont after inspecting the workspace variables' ...
      ' at different halt points.']);
disp('Opening figures can be closed.');

% detailled simulation:
clear;
% load model
disp('initializing model:');
model = convdiff_model;
disp(model);
disp('loading precomputed detailed_data:');
load('datafiles/test_rb_lin_evol_data.mat','detailed_data');
disp(detailed_data);

disp('detailed linear evolution simulation:');
model_data = gen_model_data(model);
sim_data = detailed_simulation(model, model_data);
model.plot_sim_data = @lin_evol_plot_sim_data;
params.plot_title = 'detailed simulation';
params.plot = @fv_plot;
plot_sim_data(model, model_data, sim_data, params);
disp('Please type dbcont after inspecting the workspace variables');
keyboard;

% reduced basis loading and plotting

clear('sim_data');
disp('reduced basis plotting:');
params.title = 'reduced basis vectors';
params.plot = @fv_plot;
plot_sequence(detailed_data.RB(:,1:10),detailed_data.grid,params);
%plot_mu_frequency(detailed_data.mu_values,params);
disp('Please type dbcont after inspecting the workspace variables');
keyboard;

% define complete domain for output estimation
% add some settings to perform output estimation
model.name_discfunc_element_mean = @fv_element_mean;
model.name_output_functional = 'box_mean';
model.sbox_xmin = min(detailed_data.grid.X(:));
model.sbox_xmax = max(detailed_data.grid.X(:));
model.sbox_ymin = min(detailed_data.grid.Y(:));
model.sbox_ymax = max(detailed_data.grid.Y(:));

% computation of mu and N-independent matrices 
disp('computation of mu- and N-independent offline-data:');
reduced_data = gen_reduced_data(model, detailed_data);
disp('Please type dbcont after inspecting the workspace variables');
keyboard;

% computation of mu and N-independent matrices 
disp(['Selection of N and computation of mu-independent but',...
      ' N-dependent online-data:']);
model.N = 50;
reduced_data = extract_reduced_data_subset(model, reduced_data);
disp('Please type dbcont after inspecting the workspace variables');
keyboard;

% mu-parameter variation: sets the fields of params:
%    params.c_init = 0; params.beta = 1; params.k = 5e-8; 
disp('Selection of mu and online-simulation');
model = model.set_mu(model, [0,1,5e-8]);
model.error_norm = 'l2';
disp('reduced basis simulation:');
simulation_data = rb_simulation(model, reduced_data);
disp('Please type dbcont after inspecting the workspace variables');
keyboard;

% reconstruction and plot of reduced simulation
disp('reconstruction of reduced simulation:');
simulation_data = rb_reconstruction(model,detailed_data,simulation_data);
params.title = 'reduced simulation';
params.plot = @fv_plot;
plot_sequence(simulation_data.U, detailed_data.grid, params);
figure;
rb_plot_output_estimation(simulation_data,model);
disp('Please type dbcont after inspecting the workspace variables');
keyboard;

% complete interactive gui:
clear;
disp('complete interactive gui:');
demo_rb_gui;

