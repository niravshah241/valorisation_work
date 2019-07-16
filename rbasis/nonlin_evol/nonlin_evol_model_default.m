function model=nonlin_evol_model_default

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


% localized explicit operator
model.L_E_local_ptr    = @fv_explicit_space;
% localized implicit operator
model.L_I_local_ptr    = @fv_implicit_space;

model.fv_expl_conv_weight  = 0.0;
model.fv_expl_diff_weight  = 0.0;
model.fv_expl_react_weight = 0.0;

model.fv_impl_conv_weight  = 0.0;
model.fv_impl_diff_weight  = 0.0;
model.fv_impl_react_weight = 0.0;
% implicit operators algorithm
model.implicit_operators_algorithm = @fv_operators_implicit;
model.operators_diff_implicit    = @fv_operators_zero;
model.operators_conv_implicit    = @fv_operators_zero;
model.operators_neumann_implicit = @fv_operators_zero;

% init values algorithm
model.init_values_algorithm        = @fv_init_values;
% name of function in rbmatlab/datafunc/init_values/
model.init_values_ptr = @init_values_homogeneous;
model.c_init = 0;

model.filecache_ignore_fields_in_model = {'N','Nmax','M',...
  'filecache_ignore_fields_in_model',...
  'filecache_ignore_fields_in_detailed_data',...
  'M_by_N_ratio'};
model.filecache_ignore_fields_in_detailed_data = {'RB_info'};

model.filecache_velocity_matrixfile_extract = 0;

% - assembling of decomposed initial data
model.rb_init_values                 = @rb_init_values_default;
% - detailed simulation
model.detailed_simulation            = @nonlin_evol_detailed_simulation;
% - detailed simulation with empirical interpolated discretization operator
model.detailed_ei_simulation         = @nonlin_evol_detailed_ei_simulation;
% - detailed simulation with empirical interpolated projected on rb space discretization operator
model.detailed_ei_rb_proj_simulation = @nonlin_evol_detailed_ei_rb_proj_simulation;
% - generation of model data
model.gen_model_data                 = @nonlin_evol_gen_model_data;
% - generation of detailed data
model.gen_detailed_data              = @nonlin_evol_gen_detailed_data;
% - generation of online data
model.gen_reduced_data                = @nonlin_evol_gen_reduced_data;
% - a reduced simulation
model.rb_simulation                  = @nonlin_evol_rb_simulation;
% - assembling of decomposed operators
model.rb_operators                   = @nonlin_evol_rb_operators;
% - plot detailed data
model.plot_detailed_data             = @nonlin_evol_plot_detailed_data;
model.plot_sim_data                  = @nonlin_evol_plot_sim_data;
model.rb_reconstruction              = @rb_reconstruction_default;

model.rb_init_data_basis             = @RB_init_data_basis;

model.inner_product_matrix_algorithm = @fv_inner_product_matrix;
model.inner_product                  = @fv_inner_product;
model.reduced_data_subset             = @nonlin_evol_reduced_data_subset;
model.get_inner_product_matrix       = @(detailed_data) detailed_data.W;
model.get_rb_from_detailed_data      = @(detailed_data) detailed_data.RB;
model.set_rb_in_detailed_data        = @(detailed_data,RB) ...
    setfield(detailed_data,'RB', RB);
model.get_dofs_from_sim_data         = @(sim_data) sim_data.U;
model.get_estimators_from_sim_data   = @(sim_data) sim_data.Delta;
model.get_estimator_from_sim_data    = @(sim_data) sim_data.Delta(end);
model.PCA_fixspace                   = @PCA_fixspace;
model.cached_detailed_simulation     = @cached_detailed_simulation;
model.save_detailed_simulations      = @save_detailed_simulations;
model.load_detailed_simulation       = @load_detailed_simulation;
model.set_mu                         = @set_mu_default;
model.set_time                       = @set_time_default;
model.get_mu                         = @get_mu_default;
model.get_rb_size                    = @get_rb_size;
model.laplacian_ptr                  = @(glob, U, model) U;
model.laplacian_derivative_ptr       = @(glob, U, model) ones(length(U),1);

model.l2_error_sequence_algorithm = @fv_l2_error;
model.linfty_error_sequence_algorithm = @fv_linfty_error;
model.error_algorithm = @fv_error;
model.error_norm      = 'l2l2';
model.relative_error  = false;
model.enable_error_estimator = 1;

model.RB_generation_mode     = 'greedy_uniform_fixed';
% choose extension method
model.RB_extension_algorithm = @RB_extension_PCA_fixspace;
%               'RB_extension_max_error_snapshot'};
model.RB_stop_timeout        = 2*60*60; % 2 hour

model.RB_error_indicator           = 'error';

model.RB_val_rand_seed             = 1234;
model.RB_M_val_size                = 10;
model.RB_refinement_theta          = 0.05;
model.RB_stop_max_refinement_level = 15;

model.separate_CRBs = false;

model.data_const_in_time = 1;
model.local_stencil_size = 1; % number of neighbour steps for
                              % extended eind set

model.linfty_error_sequence_algorithm = @fv_linfty_error;
model.l2_error_sequence_algorithm     = @fv_l2_error;

model.implicit_nonlinear = false;

model.newton_solver = false;
model.newton_epsilon = 1e-11;

model.geometry_transformation = 'none';

