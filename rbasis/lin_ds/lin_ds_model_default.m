function model= lin_ds_model_default
%function model= lin_ds_model_default
%
% function setting some default fields of a lin_ds model, can be
% overloaded later

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


% Bernard Haasdonk 22.9.2009

model.G_matrix_function_ptr = @(model,model_data) model_data.G;
model.set_mu = @set_mu_default;
model.get_mu = @get_mu_default;
model.set_time = @set_time_default;
model.theta = 0; % for timestepping scheme

model.detailed_simulation = @lin_ds_detailed_simulation;
model.gen_model_data = @lin_ds_gen_model_data;
model.gen_detailed_data = @lin_ds_gen_detailed_data;
model.gen_reduced_data = @lin_ds_gen_reduced_data;
model.rb_simulation = @lin_ds_rb_simulation;
model.plot_sim_data = @lin_ds_plot_sim_data;
model.plot_sim_data_state = @lin_ds_plot_sim_data_state;
model.plot_sim_data_output = @lin_ds_plot_sim_data_output;
model.plot_detailed_data = @lin_ds_plot_detailed_data;
model.reduced_data_subset = @lin_ds_reduced_data_subset;
model.rb_reconstruction = @lin_ds_rb_reconstruction;

% required for basis generation:
model.get_rb_from_detailed_data = @(detailed_data) detailed_data.V;
model.set_rb_in_detailed_data = @lin_ds_set_rb_in_detailed_data;
model.get_rb_size = @(model,detailed_data) size(detailed_data.V,2);
model.get_inner_product_matrix = ...
    @(detailed_data) detailed_data.G;
model.rb_problem_type = 'lin_ds';
model.get_estimator_from_sim_data = @(sim_data) sim_data.DeltaX(end);
model.get_dofs_from_sim_data = @(sim_data) sim_data.X;
model.get_estimators_from_sim_data = @(sim_data) sim_data.DeltaX;
model.PCA_fixspace = @PCA_fixspace;
% the following is not defined at this point...
model.init_values_algorithm = @(model,detailed_data) ...
    model.x0_function_ptr(model,detailed_data);
model.orthonormalize = @model_orthonormalize_qr;
model.inner_product = @(model,model_data,vecs1,vecs2) ...
    vecs1' * model_data.G * vecs2;
model.rb_init_data_basis = @RB_init_data_basis;

% fields to ignore, such that filecaching works in basis generation
model.filecache_ignore_fields_in_model = {'N','Nmax'}; 
model.filecache_ignore_fields_in_detailed_data = {'RB_info'};

model.enable_error_estimator = 0; % turn off as default;
model.gridtype = 'none';