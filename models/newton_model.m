function model = newton_model(params)
% function model = newton_model(params)
%
% optional fields of params:
%  model_size:   - 'big'    simulations are done on a larger grid with very
%                small time steps
%                - 'small'  simulations are done on a smaller grid with
%                bigger time steps
%  model_type:   - 'buckley_leverett' -- A Buckley-Leverett problem
%                - 'linear_diffusion' -- A linear parabolic problem with
%                    constant diffusion
%                `` \partial_{t} u - \nabla \cdot ( \mu_1 \nabla u ) = 0 ``
%                - 'exponential_diffusion' -- A non-linear parabolic problem
%                    with exponential function as diffusion coefficient
%                `` \partial_{t} u 
%                   - \nabla \cdot ( k_0 + \mu_1 u^{\mu_2} \nabla u ) ``
%                - 'eoc' -- A linear diffusion problem with known exact
%                   solution for validation of numerical scheme
%                - 'eoc_nonlinear' -- A nonlinear diffusion problem with known
%                   exact solution for validation of numerical scheme
%

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

  model = nonlin_evol_model_default;
  model.name = 'newton_model';
  model.rb_problem_type = 'newton';
  %% data that is definitely used outside of the detailed simulations
  %  model.mu_names = {'diff_m', 'diff_p', 'hill_height'};
  %  model.mu_ranges = {[0,0.5],[0.01,100],[0,0.3]};
  model.mu_names = {'diff_k0'};
  model.mu_ranges = {[0.01, 0.5]};


  %% data that might be of interest outside of the detailed simulations
  model.T  = 0.2;
  model.nt = 80;

  %% finite volume settings
%  model.diffusivity_ptr            = @diffusivity_buckley_leverett;
%  model.diffusivity_derivative_ptr = @diffusivity_buckley_leverett_derivative;
  model.diffusivity_ptr            = @diffusivity_exponential;
  model.diffusivity_derivative_ptr = @diffusivity_exponential_derivative;
  model.diff_k0                    = 0.00;
  model.diff_m                     = 1.0;
  model.diff_p                     = 2.0;

  model.num_conv_flux_ptr       = @fv_num_conv_flux_engquist_osher;
  model.num_diff_flux_ptr       = @fv_num_diff_flux_gradient;
  model.flux_linear             = 0;
  model.name_flux               = 'exponential';
  model.geometry_transformation = 'none';

  model.newton_solver     = 1;
  model.newton_steps      = 60;


  model.operators_diff_implicit    = @fv_operators_diff_implicit_gradient;
  model.operators_conv_implicit    = @fv_operators_zero;
  model.operators_neumann_implicit = @fv_operators_neumann_explicit;
  model.fv_expl_conv_weight  = 0.0;
  model.fv_impl_diff_weight  = 1.0;
  model.implicit_gradient_operators_algorithm = ...
    @fv_frechet_operators_diff_implicit_gradient;

  %    lxf_lambda               = 1.25;
  %    lxf_lambda               = fv_search_max_lxf_lambda([],    


  %% define the grid
  model = unitcube(model);
  model.xnumintervals = 100;
  model.ynumintervals = 100;

  % pointer to function in rbmatlab/datafunc/dirichlet_values
%  model.dirichlet_values_ptr = @dirichlet_values_uplow;
%  model.c_dir                = 0.2;
%  model.dir_box_xrange       = [-1 2];
%  model.dir_box_yrange       = 0.5 + [1/320 3/320];
%  model.c_dir_left           = 1;
%  model.c_dir_right          = 1;
%  model.c_dir_up             = 0;
%  model.c_dir_low            = 0.5;
%  model.dir_middle           = 0.7;
  model.dirichlet_values_ptr = @dirichlet_values_leftright;
  model.c_dir_left  = 0.1;
  model.c_dir_right = 0.6;
  model.dir_middle  = 0.5;
%  model.c_dir = 0.1;

  % pointer to function in rbmatlab/datafunc/neumann_values
  model.neumann_values_ptr = @neumann_values_homogeneous;
%  model.neumann_values_ptr = @neumann_values_outflow;
%  model.neumann_values_ptr = @neumann_values_rightflow;
  model.c_neu              = 0.00;

  % name of function in rbmatlab/datafunc/init_values/
%  model.init_values_ptr = @init_values_transformed_blobs;
  model.init_values_ptr = @init_values_bars;
%  model.init_values_ptr = @init_values_as_dirichlet;
  % parameters for data functions
  model.blob_height = 0.4;
  model.c_init      = 0.1;

  model.filecache_velocity_matrixfile_extract = 0;

  model.geometry_spline_type             = 'linear';
  model.hill_height                      = 0.0;
  model.geometry_transformation_spline_x = [ 0 0.5 1 ];
  model.geometry_transformation_spline_y = [ 0 -0.033 0 ];

  model.k       = 0.005;
  model.gravity = 0;


  % settings for CRB generation
  model.CRB_generation_mode = 'param-time-space-grid';
  model.RB_generation_mode  = 'greedy_uniform_fixed';
  model.Mmax    = 250;
  model.sMmax   = {50, 50};
  model.M       = 40;
  model.Mstrich = 5;
  model.ei_stop_on_Mmax = 1;
  % build ei-basis with WHOLE trajectory
  model.ei_target_error = 'interpol';
  model.ei_numintervals = [8];

  model.RB_stop_Nmax                = 50;
  model.RB_stop_epsilon             = 1e-9;
  model.RB_stop_max_val_train_ratio = inf;
  model.RB_numintervals            = [8];

  % target accuracy epsilon:

  % decide, whether estimator or true error is error-indicator for greedy
  model.RB_error_indicator = 'error'; % true error
  % RB_error_indicator = 'estimator'; % Delta from rb_simulation
  % RB_error_indicator = 'ei_estimator_test'; % Delta from rb_simulation testet
  % against true error

  model.divclean_mode    = false;
  model.flux_quad_degree = 1;

  model.conv_flux_ptr             = @conv_flux_linear;

  model.velocity_ptr              = @velocity_transport;
  model.transport_x               = 0.1;
  model.transport_y               = 0.1;

  model.conv_flux_derivative_ptr  = @(glob, U, params) ...
                                     params.velocity_ptr(glob, params);


  model.data_const_in_time = 0;

  % set all to dirichlet-boundary by specifying "rectangles", all
  % boundary within is set to boundary type by bnd_rect_index
  model.bnd_rect_corner1 = [-0.999999, 0; ...
                             0,       -0.9999999 ];
  model.bnd_rect_corner2 = [ 2,        0.999999; ...
                            0.999999, 2 ];
  % -1 means dirichlet, -2 means neumann
  model.bnd_rect_index   = [ -2, -2 ];

  if isfield(params, 'separate_CRBs') && params.separate_CRBs
    model.separate_CRBs = params.separate_CRBs;
    savepath_infix      = '';
  else
    savepath_infix      = '_one_CRB';
  end

  model.ei_detailed_savepath = [model.name,'_',params.model_type,...
                                savepath_infix,'_ei_data_interpol'];
  model.ei_operator_savepath = [model.name,'_',params.model_type,...
                                savepath_infix,'_ei_operators_interpol'];

  model.RB_detailed_train_savepath = model.ei_detailed_savepath;

  if isfield(params, 'model_size')
    if strcmp(params.model_size, 'big') == 1
      model.xnumintervals = 200;
      model.ynumintervals = 200;
    end
    if strcmp(params.model_size, 'small') == 1
      model.xnumintervals = 50;
      model.ynumintervals = 50;
      model.T = 2;
      model.nt = 20;
    end
  end

  if isfield(params, 'model_type')
    if strcmp(params.model_type, 'linear_diffusion')
      model.mu_names        = {'diff_k0'};
      model.mu_ranges       = {[0.01, 0.5]};
      model.ei_numintervals = [8];
      model.RB_numintervals = [8];
    elseif strcmp(params.model_type, 'buckley_leverett')
      model.bl_M                     = 0.50;
      model.bl_k                     = -1.0;
%      model.xnumintervals            = 100;
%      model.ynumintervals            = 100;
%      model.nt = 100;
      model.xnumintervals            = 50;
      model.ynumintervals            = 50;
      model.nt = 70;
      model.T                        = 0.6;
      model.conv_flux_ptr            = @conv_flux_buckley_leverett;
      model.conv_flux_derivative_ptr = @conv_flux_buckley_leverett_derivative;
%      model.num_diff_flux_ptr        = @fv_num_conv_flux_lax_friedrichs;
      model.num_conv_flux_ptr        = @fv_num_conv_flux_engquist_osher;
      model.num_diff_flux_ptr        = @fv_num_diff_flux_gradient;
      model.operators_diff_implicit  = @fv_operators_diff_implicit_gradient;
      model.implicit_gradient_operators_algorithm = ...
        @fv_operators_conv_explicit_engquist_osher;
%        @fv_operators_conv_explicit_lax_friedrichs;
      model.operators_conv_implicit    = ...
        @fv_operators_conv_explicit_engquist_osher;
%      model.operators_conv_explicit    = ...
%        @fv_operators_conv_explicit_engquist_osher;
%        @fv_operators_conv_explicit_lax_friedrichs;
      model.fv_impl_conv_weight      = 1.0;
%      model.fv_expl_conv_weight      = 1.0;
      model.fv_impl_diff_weight      = 1.0;
      model.neumann_values_ptr = @neumann_values_rightflow;
      model.lxf_lambda               = 0.3;
      model.init_values_ptr          = @init_values;
      model.blob_height = 0.5;
      model.init_values_ptr          = @init_values_as_dirichlet;
      model.dirichlet_values_ptr     = @init_values;
      model.diffusivity_ptr = @diffusivity_homogeneous;
      model.k = 0.0000;
%      model.c_dir = -0.1;
      model.c_init = 0.1;
      model.c_dir_left = 0.2;
      model.c_dir_right = 0.2;
      model.yrange = [0.2 1];
      model.ei_time_splits = 5;
      model.Mmax           = 300;
      model.RB_stop_timeout = inf;
%      model.dirichlet_values_ptr     = @dirichlet_values_quarter_of_five;
%      model.c_dir_right = 0.0;
%      model.c_dir_left  = 0.6;
%      model.bnd_rect_corner1 = [-0.1, -0.1; ...
%                                 0.9, 0.9; ...
%                                -0.1, 0.099; ...
%                                 0.099, -0.1 ]';
%      model.bnd_rect_corner2 = [ 0.1,        0.1; ...
%                                 1.1,        1.1; ...
%                                 0.901,      1.1; ...
%                                 1.1,        0.901 ]';
%      model.bnd_rect_index           = [-1, -1, -2, -2];
      model.conv_a = 1.0;
      model.bnd_rect_index   = [-1, -1];
%      model.mu_names = {'k','blob_height','c_init','conv_a','bl_M'};
%      model.mu_ranges = { [0 1e-3], [0.1 0.5], [-0.1 0.1], [0 1], [0.2, 0.5] };
      model.mu_names = {'blob_height','conv_a','bl_k'};
      model.mu_ranges = { [0.2 0.5], [0.2 0.5], [-1, -0.5]};
      model.ei_numintervals = [6, 6, 6];
      model.RB_numintervals = [6, 6, 6];
      % purely implicit scheme: only interpolate the implicit operator therefore!
      model.ei_space_operators = { model.L_I_local_ptr };
    elseif strcmp(params.model_type, 'exponential_diffusion')
      model.ei_numintervals = [6,10];
      model.RB_numintervals = [6,10];
      model.diff_k0         = 0.0000;
      model.diff_m          = 0.01;
      model.c_init = 0.0001;
      model.newton_regularisation = 0;
      model.diff_p          = 2.0;
      model.T = 1;
      model.mu_names        = {'diff_m','diff_p','diff_k0'};
      model.mu_ranges       = {[0.0,0.01],[1,5],[0,0.0001]};
      model.ei_numintervals = [4,6,3];
      model.RB_numintervals = [4,6,3];
      model.nt              = 50;
      model.RB_stop_Nmax    = 50;
      model.RB_stop_epsilon = 1e-6;
      model.RB_stop_timeout = inf;
      model.bnd_rect_index  = [-2,-2];
      model.separate_CRBs   = false;
      model.ei_space_operators = { model.L_I_local_ptr };

      if isfield(params, 'use_laplacian') && params.use_laplacian
        model.laplacian_derivative_ptr = @(glob, U, model) ...
          real((model.diff_m*U.^(model.diff_p-1)));
        model.laplacian_ptr            = @(glob, U, model) ...
          real((1/model.diff_p * model.diff_m * U.^(model.diff_p)));
        model.diffusivity_ptr          = @(glob, U, model) ...
          struct('K',1,'epsilon',0);
        model.diffusivity_derivative   = @(glob, U, model) ...
          struct('K',1,'epsilon',0);
        model.implicit_gradient_operators_algorithm = ...
          @fv_operators_diff_implicit_gradient;
      end
    elseif strcmp(params.model_type, 'eoc')
      model.xnumintervals        = 10;
      model.ynumintervals        = 10;
      model.diff_k0              = 0.05;
      model.diff_m               = 0;
      model.nt                   = 50;
      model.init_values_ptr      = @exact_function_heat_equation;
      model.dirichlet_values_ptr = @exact_function_heat_equation;
      model.bnd_rect_index       = [ -1, -1 ];
    elseif strcmp(params.model_type, 'eoc_nonlinear')
      model.xnumintervals        = 10;
      model.ynumintervals        = 10;
      model.newton_regularisation = 0;
      model.diff_k0              = 0.0;
      model.diff_m               = 1;
      model.diff_p               = 1.0;
      model.nt                   = 150;
      model.T                    = 1.0;
      model.init_values_ptr      = @exact_function_plaplace;
      model.dirichlet_values_ptr = @exact_function_plaplace;
      model.data_const_in_time   = false;
      model.bnd_rect_index       = [ -1, -2 ];
    elseif strcmp(params.model_type, 'buckley_leverett')
      error('buckley leverett problem type is not implemented yet.');
    end
  end

  model.model_type = params.model_type;

  model.implicit_nonlinear = true;

  model.data_const_in_time = false;

  model = model_default(model);

  if ~isfield(model, 'ei_space_operators')
    model.ei_space_operators = { model.L_E_local_ptr, model.L_I_local_ptr };
  end

  if isfield(params, 'verbose')
    model.verbose = params.verbose;
  else
    model.verbose = 0;
  end
  model.debug   = 0;

end

function U0=init_values(glob, params)
  if ~isfield(params, 'blob_radius')
    params.blob_radius = 0.2;
  end
  params.blob_height = params.blob_height - params.c_init;
  if ~isfield(params, 'blob_height2')
    params.blob_height2 = params.blob_height;
  end
  params.blob_height2 = params.blob_height2 - params.blob_height;

  decomp_mode = params.decomp_mode;
  if decomp_mode == 2
    U0 = [params.c_init params.blob_height params.blob_height2];
  else
    X = glob(:,1);
    Y = glob(:,2);
    B1 = sqrt((X(:)-0.75).^2+(Y(:)-0.75).^2);
    if decomp_mode == 0
      U0 = params.c_init * ones(size(X));
      U0 = U0 + params.blob_height * (B1 < params.blob_radius);
      U0 = U0 + params.blob_height2 * (B1 < 1/2*params.blob_radius);
    elseif decomp_mode == 1
      U0 = cell(3,1);
      U0{1} = ones(size(X(:)));
      U0{2} = (B1 < params.blob_radius);
      U0{3} = (B1 < 1/2*params.blob_radius);
    end
  end
end

% vim: set et sw=2:
%| \docupdate 
