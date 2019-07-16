function model = riemann_burgers_model(params)
%function model = riemann_burgers_model(params)
%
% function generating a model of a non-viscous burgers flow
% developing either a rarefaction wave or a moving shock-front.
% Model is used for the hyp2008 paper and demo_rb_riemann_burgers
% Currently, params is ignored.

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


% Bernard Haasdonk 18.3.2010

if nargin < 1
  params = [];
%  params.size = 1;
end

model = nonlin_evol_model_default;

model.gridtype = 'rectgrid';
if isfield(params,'xnumintervals')
  model.xnumintervals = params.xnumintervals;
else
  model.xnumintervals = 100;
end;

if model.xnumintervals>100
  error('please adjust time-step number to fine grid resolution!');
end;

if isfield(params,'xnumintervals')
  model.ynumintervals = params.ynumintervals;
else
  model.ynumintervals = 100;
end;

model.xrange = [0,1];
model.yrange = [0,1];

% set upper and lower to neuman noflow
% set left and right to dirichlet

%epsilon = 1e-6;
model.bnd_rect_corner1 = [-1, 0-eps, 1-eps; ...
		           -1, 0-eps, 0-eps];
model.bnd_rect_corner2 = [ 2, 0+eps, 1+eps; ...
		            2, 1+eps, 1+eps];

% -1 means dirichlet, -2 means neuman
model.bnd_rect_index =   [-2,-1,-1]; 

% time discretization
%model.T = 0.5/10;
%model.nt = 100/10; % guess and adjust later
model.T = 0.5;
model.nt = 100; % guess and adjust later

%model.T = 0.4;
%model.nt = 80; % guess and adjust later
%model.nt = 100; % guess and adjust later
%model.T = 0.1;
%model.nt = 50; % guess and adjust later
%model.400 => non-bounded u

% output settings
model.verbose = 9;
model.axis_equal = 1;
model.clim = [-1,1];
model.debug = 0;

% data functions
%model.name_init_values = 'blobs'; 
%model.radius = 0.1;
%model.gamma = 100;
%model.beta = 1; % init data weight between 0 and 1

% name of function in rbmatlab/datafunc/init_values
%model.name_init_values = 'waveproduct'; 
%model.name_init_values = 'leftright'; 
%model.name_init_values = 'as_dirichlet'; 
model.init_values_ptr = @init_values_as_dirichlet; 

%model.c_init_max = 1.0;
%model.c_init_phase_x = 0.0;
%model.c_init_phase_y = 0.0;
%model.c_init_freq_x = 2*pi;
%model.c_init_freq_y = 2*pi;

%model.c_init_lo = 1.0;
%model.c_init_lo = 0.0;
%model.c_init_lo = -1.0;

% name of function in RBmatlab/datafunc/dirichlet_values
%disp('todo: dirichlet values!');
%model.name_dirichlet_values = 'leftright';
model.dirichlet_values_ptr = @dirichlet_values_leftright;
model.c_dir_left = -1.0;
model.c_dir_right = 0.0;
model.dir_middle = 0.5;

% name of function in RBmatlab/datafunc/neuman_values
model.neumann_values_ptr = @neumann_values_zero;
%model.name_neuman_values = 'zero';
%model.name_neuman_values = 'homogeneous';
%model.c_neu = 0;

% convective flux
%model.name_flux = 'burgers_parabola';
model.conv_flux_ptr = @conv_flux_burgers;
%| \todo this could also be done explicitly (analytic derivative)
model.conv_flux_derivative_ptr = @conv_flux_forward_difference;
% for simplifying computation in case of linear flux:
model.flux_linear = 0;
% diagonal velocity field
model.flux_vx = 1;
model.flux_vy = 0;
%model.flux_quad_degree = 2;
model.flux_pdeg = 2; % burgers exponent

% rb-formulation
%model.mu_names = {'c_init_lo','vrot_angle'}; 
% themiddle parameter makes u_init not decomposable...
%model.mu_names = {'c_dir_left','c_dir_right','dir_middle'}; 
model.mu_names = {'c_dir_left','c_dir_right','flux_vx'}; 
%model.mu_ranges = {[0,1],[0,1],[0,1]};
model.mu_ranges = {[-1,1],[-1,1],[-1,1]};
model.rb_problem_type = 'nonlin_evol';
%model.problem_type = 'lin_evol';

% finite volume settings
%model.name_diffusive_num_flux = 'none';
%model.name_convective_num_flux = 'lax-friedrichs';
%model.lxf_lambda = 1.66
%model.lxf_lambda = 1.25;
%model.lxf_lambda = 0.35355;
%model.lxf_lambda = fv_search_max_lxf_lambda([],model);
%keyboard;
%model.name_convective_num_flux = 'engquist-osher';
model.num_conv_flux_ptr = @fv_num_conv_flux_engquist_osher;
model.fv_expl_conv_weight = 1.0;

% projection of analytical initial data to discrete function
model.init_values_algorithm = @fv_init_values;
%model.implicit_operators_algorithm = @fv_operators_implicit;
% get mass matrix for inner product computation from DOF-vectors
model.data_const_in_time = 1;

%model.init_values_algorithm = 'fv_init_values';
%model.implicit_operators_algorithm = 'fv_operators_implicit';
% get mass matrix for inner product computation from DOF-vectors
%model.inner_product_matrix_algorithm = 'fv_inner_product_matrix';
%model.data_const_in_time = 1;

% settings for empirical interpolation
par.range = model.mu_ranges;
model.ei_numintervals = [4,4,4]; % 5^3 parameter grid
%model.ei_numintervals = [9]; % 8x8 parameter grid
%model.Mmax = 300;
model.Mmax = 100;
model.M = 100;
model.ei_detailed_savepath = 'nonlin_symmetry_detailed_200';
model.ei_operator_savepath = 'nonlin_symmetry_ei_operator_200'; 
model.ei_time_indices = 1:model.nt; 
model.CRB_generation_mode = 'param-time-space-grid';
model.ei_target_error = 'interpol';
model.flux_quad_degree = 1;


% the following may be required:
%model.L_E_local_name = 'fv_conv_explicit_space';

%model.ei_numintervals = [2,2,2]; % 5^3 parameter grid
%model.ei_numintervals = [4,4]; % 5^2 parameter grid
%model.Mmax = 300;
%model.Mmax = 150;
%model.ei_detailed_savepath = ei_detailed_savepath;
%model.ei_operator_savepath = ei_operator_savepath;
%model.ei_time_indices = 1:model.nt; 
%model.CRB_generation_mode = 'param-time-space-grid';
%model.ei_target_error = 'interpol';

% settings for reduced basis generation
%model.detailed_simulation_algorithm = 'detailed_simulation';
%model.operators_algorithm = 'fv_operators_implicit_explicit';
%model.init_values_algorithm = 'fv_init_values';
%model.inner_product_matrix_algorithm = 'fv_inner_product_matrix';
%model.error_algorithm = 'fv_error';
%model.lxf_lambda = 1.0194e+003;
%model.data_const_in_time = 1;
model.error_norm = 'l2';
%model.RB_extension_algorithm = 'RB_extension_PCA_fixspace';
% 		    'RB_extension_max_error_snapshot'}; 
%model.RB_stop_timeout = 2*60*60; % 2 hours		
%model.RB_stop_epsilon = 1e-5; 
%model.RB_error_indicator = 'error'; % true error
%model.RB_error_indicator = 'estimator'; % Delta from rb_simulation
%model.RB_stop_Nmax = 100; 
%model.RB_generation_mode = 'uniform_fixed';
model.RB_generation_mode = 'PCA_trajectories';
model.RB_mu_list = {[1,0,1],[0,1,-1]}; % two shock-trajectories in
                                        % different direction
%model.RB_numintervals = model.ei_numintervals;
%model.RB_detailed_train_savepath = model.ei_detailed_savepath;

% parameters such that demo_rb_gui looks nicer:
model.yscale_uicontrols = 0.6;
model.xscale_gui = 0.5;
model.show_colorbar = 0;
model.clim = [-1,1];

model.t               = 0;
%model.tstep           = 1;
model.decomp_mode     = 0;
model.orthonormalize  = @model_orthonormalize_qr;
model.plot            = @fv_plot;
model.error_estimation = 0; % turn off rb error estimation

%model.gen_reduced_data = @nonlin_evol_gen_reduced_data;
%model.rb_operators = @nonlin_evol_rb_operators;
%model.set_time = @set_time_default;
%model.set_mu = @set_mu_default;
%model.reduced_data_subset = @nonlin_evol_reduced_data_subset;
%model.rb_simulation = @nonlin_evol_rb_simulation;
