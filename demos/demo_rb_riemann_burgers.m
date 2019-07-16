%demo_rb_riemann_burgers
%
% Script demonstrating the RB-model for the discontinuous 
% initial data for the nonlinear burgers equation. 
% Both moving shock front and rarefaction waves can be modelled.
% discretization is with explicit fv scheme and empirical
% interpolation of the space operator.
% The model demonstrates the space-dimension reduction: the 
% subgrid that is extracted in the offline phase is exactly the row
% of cells along the bottom line.
%
% To start: set u_left or u_right to some different value.
% Then observe either a moving shock or a rarefaction wave
% depending on the direction of the velocity.

% Bernard Haasdonk 10.9.2008

help demo_rb_riemann_burgers

detailedfname = 'riemann_burgers_detailed.mat';
load(detailedfname);
%model = riemann_burgers_model;
model.newton_solver = false;
model.implicit_nonlinear = false;

detailed_data.grid = construct_grid(model);

%[model,detailed_data,plot_params] = renew_model(model,detailed_data);

%model = riemann_burgers_model;

%model.N = model.get_rb_size(model,detailed_data);
%size(detailed_data.RB,2);
%model.M = size(detailed_data.QM{1},2);

rb_plot_interpolation_points(detailed_data,model);

%model.gen_reduced_data = @nonlin_evol_gen_reduced_data;
%model.rb_operators = @nonlin_evol_rb_operators;
%rb_plot_interpolation_points(detailed_data,model);
%model.set_time = @set_time_default;
%model.set_mu = @set_mu_default;
%model.reduced_data_subset = @nonlin_evol_reduced_data_subset;
%%model.rb_simulation = @nonlin_evol_rb_simulation;

% temporary: transform to "model"
% [model, detailed_data, plot_params] = renew_model(params, detailed_data);

demo_rb_gui(model,detailed_data); %,[],plot_params,'riemann burgers');
