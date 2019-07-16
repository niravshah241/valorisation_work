function model = convdiff_dune_model(params)

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


model.model_name = 'convdiff_dune_model';
model.rb_problem_type = 'lin_evol';
%% data that is definitely used outside of the detailed simulations
%% These values are now set in dunerbconvdiff( 'init_model')

model.L_I_inv_norm_bound = 1;
model.L_E_norm_bound     = 1;

model.error_algorithm = @fv_error;
model.error_norm      = 'l2';
model.relative_error  = false;

%% reduced basis generation settings

% choose generation mode
model.RB_generation_mode = 'greedy_uniform_fixed';
% alternative:
% RB_generation_mode = 'random_fixed';

% choose extension method
model.RB_extension_algorithm = @dune_RB_extension_PCA_fixspace;
% stop timeout
model.RB_stop_timeout = 60*60; % 1 hour
% stop on maximum number of base functions
model.RB_stop_Nmax = 30;
% stop on epsilon
model.RB_stop_epsilon = 1e-9; % 0.003 is already realized by init data!!
% stop on overfitting (needs a validation test set)
model.RB_stop_max_val_train_ratio = inf;

%model.RB_test_size = 1000;

% fix the rand seed in case of RB_generation_mode = random_fixed
model.RB_train_rand_seed = 0.5;
% trainig set size in case of RB_generation_mode = random_fixed
model.RB_train_size = 50;

% choose the error indicator
model.RB_error_indicator = 'estimator'; % Delta from rb_simulation
% alternative:
% model.RB_error_indicator = 'error'; % true error

model.xrange        = [0,1000e-6]; % grid x-coordinate range
model.yrange        = [0,200e-6];  % grid y-coordinate range
model.xnumintervals = 200;  % number of grid cellls in x-direction
model.ynumintervals = 40;   % number of grid cellls in y-direction

model.t = 0;
model.decomp_mode = 0;

model.coeff_ops = [];

model.data_const_in_time = 0;

model.detailed_simulation = @dune_detailed_simulation;
model.gen_detailed_data   = @lin_evol_gen_detailed_data;
model.gen_reduced_data     = @lin_evol_gen_reduced_data;
model.rb_init_data_basis  = @dune_rb_init_data_basis;
model.gen_model_data      = @dune_gen_model_data;
model.rb_init_values      = @dune_rb_init_values;
model.rb_operators        = @dune_lin_evol_rb_operators;
model.rb_simulation       = @lin_evol_rb_simulation;
model.rb_reconstruction   = @dune_lin_evol_rb_reconstruction;
model.set_mu              = @dune_set_mu;
model.set_time            = @set_time;
model.get_rb_size         = @get_rb_size;
model.get_rb_from_detailed_data = @(detailed_data) detailed_data.RB;
model.set_rb_in_detailed_data = @(detailed_data, RB) ...
    setfield(detailed_data, 'RB', RB);
model.get_dofs_from_sim_data    = @(sim_data) sim_data.U;
model.reduced_data_subset        = @lin_evol_reduced_data_subset;
model.get_estimator_from_sim_data = @(sim_data) sim_data.Delta(end);
model.get_mu    = @dune_get_mu;

%orthonormalize = @orthonormalize_gram_schmidt;
%inner_product  = @fv_inner_product;
%cached_detailed_simulation = @cached_detailed_simulation;
%save_detailed_simulations = @save_detailed_simulations;
%load_detailed_simulation = @load_detailed_simulation;

%mexptr = @dunerbconvdiff;
model.mexptr = @mexclient;

model.verbose = 1;
model.debug = 0;

model.data_const_in_time = 1;

% server_params.serverhost = 'rbevol.uni-muenster.de'
% server_params.port       = '1909'
% mexclient('init_server', server_params)

addpath( fullfile( getenv('DUNERBHOME'), 'mexclient' ) );
addpath( getenv('DUNERBHOME') );

tmp_struct = model.mexptr('init_model');

%model.diffusion = 0;
%model.velocity

model.T         = tmp_struct.T;
model.nt        = tmp_struct.nt;
model.dt        = model.T / model.nt;
model.mu_names  = tmp_struct.mu_names;
model.mu_ranges = tmp_struct.mu_ranges;
model.mu        = model.mexptr('get_mu');


tmp = model.mexptr('rb_symbolic_coefficients');
sym_coeffs = cell(1,4);
for i = 1:4
  sym_coeffs{i} = '@(mu) [ ';
  for j = 1:size(tmp{i},2);
    if j > 1
      sym_coeffs{i} = [ sym_coeffs{i}, ', ' ];
    end
    sym_coeffs{i} = [ sym_coeffs{i}, tmp{i}{j} ];
  end
  sym_coeffs{i} = [ sym_coeffs{i} , ' ]' ];
end

sym_coeffs = cellfun(@eval, sym_coeffs, 'UniformOutput', false);

model.coeff_ops.LL_I_ptr = sym_coeffs{1};
model.coeff_ops.LL_E_ptr = sym_coeffs{2};
model.coeff_ops.bb_ptr   = sym_coeffs{3};
model.coeff_ops.u0_ptr   = sym_coeffs{4};

model.RB_numintervals = 5 * ones(size(model.mu_ranges));

function [mu] = dune_get_mu(model)
  mu = model.mexptr('get_mu');
%  if any(mu(:) ~= model.mu(:))
%    warning(sprintf('mu values in RBmatlab and dune-rb differ! [%s != %s]', mat2str(mu(:)), mat2str(model.mu(:)) ) );
%  end

function reconstruct_and_compare(model, rb_sim_data)
  model.mexptr('reconstruct_and_compare', rb_sim_data.a);
%| \docupdate 
