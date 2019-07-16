function model=lin_evol_model_default

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


model.rb_problem_type = 'lin_evol';

model.L_I_inv_norm_bound = 1;
model.L_E_norm_bound     = 1;

model.operators_ptr               = @fv_operators_implicit_explicit;
model.operators_diff_implicit     = @fv_operators_zero;
model.operators_diff_explicit     = @fv_operators_zero;
model.operators_conv_implicit     = @fv_operators_zero;
model.operators_conv_explicit     = @fv_operators_zero;
model.operators_neumann_explicit  = @fv_operators_zero;
model.operators_neumann_implicit  = @fv_operators_zero;

model.detailed_simulation         = @lin_evol_detailed_simulation;
model.gen_detailed_data           = @lin_evol_gen_detailed_data;
model.gen_reduced_data            = @lin_evol_gen_reduced_data;
model.rb_init_data_basis          = @RB_init_data_basis;
model.gen_model_data              = @lin_evol_gen_model_data;
model.rb_init_values              = @rb_init_values_default;
model.rb_operators                = @lin_evol_rb_operators;
model.rb_simulation               = @lin_evol_rb_simulation;
model.rb_reconstruction           = @rb_reconstruction_default;
model.plot                        = @fv_plot;
model.plot_detailed_data          = @lin_evol_plot_detailed_data;
model.plot_sim_data               = @lin_evol_plot_sim_data;
model.get_rb_size                 = @get_rb_size;
model.set_mu                      = @set_mu_default;
model.set_time                    = @set_time_default;
model.get_mu                      = @get_mu_default;
model.get_inner_product_matrix    = @(detailed_data) detailed_data.W;
model.get_rb_from_detailed_data   = @(detailed_data) detailed_data.RB;
model.set_rb_in_detailed_data     = @(detailed_data,RB) ...
                                    setfield(detailed_data, 'RB', RB);

model.get_estimator_from_sim_data = @(sim_data) sim_data.Delta(end);
model.get_dofs_from_sim_data      = @(sim_data) sim_data.U;
model.get_estimators_from_sim_data = @(sim_data) sim_data.Delta;
model.PCA_fixspace                = @PCA_fixspace;
model.orthonormalize              = @model_orthonormalize_gram_schmidt;
model.inner_product               = @(model,model_data,vecs1,vecs2) ...
                                      vecs1' * model_data.W * vecs2;
model.reduced_data_subset         = @lin_evol_reduced_data_subset;

model.mass_matrix                 = @fv_mass_matrix;
model.divclean_mode               = 0;

model.l2_error_sequence_algorithm = @fv_l2_error;
model.linfty_error_sequence_algorithm = @fv_linfty_error;
model.error_algorithm = @fv_error;
model.error_norm      = 'l2';
model.error_estimation = 1; % turn on rb error estimation as default

model.affinely_decomposed = true;

model.compute_output_functional = 0;

model.data_const_in_time = 1;

%    model.lxf_lambda = 0.1172; % lambda value for Lax-Friedrichs diffusivity
%    model.lxf_lambda = 0.001172; % lambda value for Lax-Friedrichs diffusivity
model.lxf_lambda = 1.0194e+003;

model.enable_error_estimator = 1;
