function model = lin_stat_model_default
%function model = lin_stat_model_default
%
% function initializing some fields of a lin-stat model

% This program is open source.  For license terms, see the COPYING file.
%
% --------------------------------------------------------------------
% ATTRIBUTION NOTICE:
% This product includes software developed for the RBmatlab project at
% (C) Universities of Stuttgart and Münster, Germany.
%
% RBmatlab is a MATLAB software package for model reduction with an
% emphasis on Reduced Basis Methods. The project is maintained by
% M. Dihlmann, M. Drohmann, B. Haasdonk, M. Ohlberger and M. Schaefer.
% For Online Documentation and Download we refer to www.morepas.org.
% --------------------------------------------------------------------


% B. Haasdonk 21.2.2011

model = [];
model.verbose = 5;
model.debug = 0;
model.rb_problem_type = 'lin_stat';

model.mu_names = {};
model.mu_ranges = {};
model.set_mu = @set_mu_default;
model.get_mu = @get_mu_default;
model.decomp_mode = 0; % default: complete evaluation
model.gen_model_data = @lin_stat_gen_model_data;
model.detailed_simulation = @lin_stat_detailed_simulation;
model.plot_sim_data = @lin_stat_plot_sim_data;

model.gen_detailed_data = @lin_stat_gen_detailed_data;
model.gen_reduced_data = @lin_stat_gen_reduced_data;
model.reduced_data_subset = @lin_stat_reduced_data_subset;
model.rb_simulation = @lin_stat_rb_simulation;
model.rb_reconstruction = @lin_stat_rb_reconstruction;
model.compute_output_functional = 0;
model.operators = @fem_operators;
model.operators_output = @fem_operators_output;

model.get_dofs_from_sim_data = @(sim_data) sim_data.uh.dofs;
model.get_inner_product_matrix = @(detailed_data) ...
    detailed_data.df_info.regularized_h10_inner_product_matrix;
model.RB_generation_mode = 'lagrangian';
model.set_rb_in_detailed_data=@lin_stat_set_rb_in_detailed_data;
model.get_rb_size= @(model,detailed_data) size(detailed_data.RB,2);
model.get_estimators_from_sim_data= @(sim_data) sim_data.Delta;
% for demo_rb_gui:
model.is_stationary = 1;
model.axes_tight = 1;

model.has_diffusivity = 0;
model.has_source = 0;
model.has_reaction = 0;
model.has_advection = 0;
model.has_output_functional = 0;
model.has_dirichlet_values = 0;
model.has_neumann_values = 0;
model.has_robin_values = 0;

