function model = richards_fv_model(params)
% function model = richards_fv_model(params)
%
% optional fields of params:
%  model_size:   - 'big'    simulations are done on a larger grid with very
%                small time steps
%                - 'small'  simulations are done on a smaller grid with
%                bigger time steps
%  model_type:   
%                - 'linear_heat_trapezoidal'   on transformed domain linear heat equation
%                with a trapezoidal geometry parametrisation. on reference
%                domain a linear convection-diffusion-reaction problem. It
%                still depends on empirical interpolation, however, because
%                of non-affine parameter dependencies of the operators
%                introduced by geometry parametrisation. Note, that
%                diffusion term is discretized explicitly.
%                - 'linear_heat_polynomial'  on transformed domain linear heat equation
%                with polynomial geometry parametrization. On reference
%                domain a linear convection-diffusion-reaction problem.
%                Note, that diffusion term is discretized explicitly.
%                - 'nonlinear_richards_trapezoidal'  on transformed domain real Richards
%                equation example with trapezoidal geometry parametrisation
%                and a gravity induced convection term. On reference
%                domain a non-linear convection-diffusion-reaction problem.
%                Note, that diffusion term is discretized explicitly.
%                - 'linear_heat_polynomial_implicit'  On transformed domain
%                lienar heat equation with polynomial geometry
%                parametrisation. On reference domain a linear
%                convection-diffusion-reaction problem with non-affine
%                parameter dependence (resolved by EI). The diffusion term
%                is discretized implicitly. 
%                - 'test'  very small heat equation example suitable for
%                quick test runs

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
  model.name = 'richards_fv_model';
  model.rb_problem_type = 'nonlin_evol';
  %% data that is definitely used outside of the detailed simulations
  %    mu_names = {'c_init','hill_height'};
  %    mu_ranges = {[-1,0],[-0.1,0.1]};
  model.mu_names = {'c_dir_low','hill_height'};
  model.mu_ranges = {[0,0.5],[0.0,0.4]};


  %% data that might be of interest outside of the detailed simulations
  model.T  = 0.9;
  model.nt = 300;

  %% finite volume settings
  model.diffusivity_ptr        = @diffusivity_homogeneous;
  model.diffusivity_tensor_ptr = @diffusivity_tensor_richards;
  %    diffusivity_ptr = @diffusivity_linear_gradient;
  %    diff_left  = 0.001;
  %    diff_right = 0.002;


  model.num_conv_flux_ptr = @fv_num_conv_flux_engquist_osher;
  model.num_diff_flux_ptr = @fv_num_diff_flux_gradient_tensor;
  model.flux_linear       = true;
  model.name_flux         = 'richards';


  model.operators_diff_implicit    = @fv_operators_diff_implicit_gradient_tensor;
  model.operators_conv_implicit    = @fv_operators_zero;
  model.operators_neumann_implicit = @fv_operators_zero;
  model.fv_expl_conv_weight  = 1.0;
  model.fv_expl_diff_weight  = 1.0;
  model.fv_expl_react_weight = 1.0;

  %    lxf_lambda               = 1.25;
  %    lxf_lambda               = fv_search_max_lxf_lambda([],    


  %% define the grid
  model = unitcube(model);
  model.xnumintervals = 100;
  model.ynumintervals = 80;

  % pointer to function in rbmatlab/datafunc/dirichlet_values
  model.dirichlet_values_ptr = @dirichlet_values_uplow;
  model.c_dir                = 0.1;
  model.dir_box_xrange       = [-1 2];
  model.dir_box_yrange       = 0.5 + [1/320 3/320];
  model.c_dir_left           = 1;
  model.c_dir_right          = 1;
  model.c_dir_up             = 0;
  model.c_dir_low            = 0.5;
  model.dir_middle           = 0.7;
  %    name_dirichlet_values = 'homogeneous';
  %    c_dir = 0;

  % pointer to function in rbmatlab/datafunc/neumann_values
  %neumann_values_ptr = @neumann_values_homogeneous;
  model.neumann_values_ptr = @neumann_values_zero;
  model.c_neu             = 0;

  % name of function in rbmatlab/datafunc/init_values/
  model.init_values_ptr = @init_values_transformed_blobs;
  % parameters for data functions
  model.blob_height = 0.1;

  model.filecache_velocity_matrixfile_extract = 0;

  model.geometry_spline_type             = 'cubic';
  model.hill_height                      = 0.0;
  model.geometry_transformation_spline_x = [ 0 0.5 1 ];
  model.geometry_transformation_spline_y = [ 0 -0.033 0 ];

  model.k       = 0.005;
  model.gravity = 0;


  % settings for CRB generation
  model.CRB_generation_mode = 'param-time-space-grid';
  model.Mmax    = 150;
  model.sMmax   = {50, 50};
  model.M       = 40;
  model.Mstrich = 5;
  model.ei_stop_on_Mmax = 1;
  % build ei-basis with WHOLE trajectory
  model.ei_target_error = 'interpol';
  model.ei_numintervals = [2,4];

  model.RB_stop_Nmax                = 50;
  model.RB_stop_epsilon             = 1e-5;
  model.RB_stop_max_val_train_ratio = inf;
  model.RB_numintervals            = [2,4];

  % target accuracy epsilon:

  % decide, whether estimator or true error is error-indicator for greedy
  model.RB_error_indicator = 'error'; % true error
  % RB_error_indicator = 'estimator'; % Delta from rb_simulation
  % RB_error_indicator = 'ei_estimator_test'; % Delta from rb_simulation testet against true error

  model.divclean_mode    = false;
  model.flux_quad_degree = 1;

  model.conv_flux_ptr             = @conv_flux_linear;

  model.velocity_ptr              = @velocity_richards;

  model.conv_flux_derivative_ptr  = @(glob, U, params) params.velocity_ptr(glob, params);


  model.data_const_in_time = 0;

  % set all to dirichlet-boundary by specifying "rectangles", all
  % boundary within is set to boundary type by bnd_rect_index
  model.bnd_rect_corner1 = [-0.999999, 0; ...
                             0,       -0.9999999 ];
  model.bnd_rect_corner2 = [ 2,        0.999999; ...
                            0.999999, 2 ];
  % -1 means dirichlet, -2 means neumann
  model.bnd_rect_index   = [ -1, -2 ];

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

  model.ei_space_operators      = { model.L_E_local_ptr };

  model.geometry_transformation  = 'spline';

  model.stencil_mode        = 'vertex';
  model.local_stencil_size  = 2;

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
    if strcmp(params.model_type, 'nonlinear') == 1
      model.geometry_spline_type = 'cubic';
      model.nt = 300;
    elseif strcmp(params.model_type, 'linear_heat_trapezoidal') == 1
      model.geometry_spline_type = 'affine';
      model.T  = 0.5;
      model.nt = 260;
      model.c_init = 0;
    elseif strcmp(params.model_type, 'richards_affine') == 1
      model.geometry_spline_type = 'affine';
      model.diffusivity_ptr      = @diffusivity_richards_nonlinear;
      model.init_values_ptr      = @init_values_transformed_blobs_richards;
      % set parameters for nonlinear gradient
      model.k                    = 0.002;
      model.gravity              = 0.5;
      model.clim                 = [0.22,0.40];
      % set new time parameters
      model.T                    = 0.01;
      model.nt                   = 195;
      model.c_dir_up             = 0.32;
      model.c_dir_low            = 0.35;
      model.xnumintervals        = 80;
      model.ynumintervals        = 100;
      % set default mu values
      model.c_init               = 0.34;
      model.hill_height          = 0.50;
      % choose different mu vector
      model.mu_names             = {'hill_height','c_init'};
      model.mu_ranges            = {[0,0.5],[0.22,0.34]};
      % only neumann boundaries
      model.bnd_rect_index  = [ -2, -2 ];

      model.ei_numintervals        = [5,6];
      model.RB_numintervals        = [5,6];
      ht                           = @(t)(t-0.218)./(0.52-0.218);
      model.richards_perm_ptr      = @(t)(2.1*ht(t).^0.5.*(1-(1-ht(t).^(7/6)).^(6/7)).^2);
      model.richards_retention_ptr = @(t)(35.2574./(1./(3.311258.*t-0.7218543).^(6/7)-1).^(6/7)./(3.311258.*t-0.721854).^(13/7));
      model.postprocess            = @postprocess_gravity;
      model.RB_stop_Nmax           = 10;
    elseif strcmp(params.model_type, 'implicit_nonaffine_linear') == 1
      model.implicit_nonlinear   = true;
      model.hill_height          = 0.4;
      model.geometry_spline_type = 'cubic';
      %            model.mu_names             = {'hill_height', 'blob_height'};
      model.mu_ranges            = {[0,0.4], [0.0,0.3]};
      %            model.c_init               = 0;
      %            model.bnd_rect_index  = [ -2, -2 ];
      model.k                    = 0.005;
      %    k                     = 0.005;
      model.nt                   = 100;
      model.fv_expl_diff_weight  = 0.0;
      model.fv_impl_diff_weight  = 1.0;
      model.data_const_in_time   = 1;
      model.T                    = 1.5;
      model.ei_space_operators      = { model.L_E_local_ptr, model.L_I_local_ptr };
    elseif strcmp(params.model_type, 'test')
      model.T                    = 0.15;
      model.nt                   = 30;
      model.Mmax                 = 20;
      model.RB_stop_Nmax         = 5;
      model.M                    = model.Mmax;
      model.Mstrich              = 0;
      model.geometry_spline_type = 'affine';
      model.data_const_in_time   = true;
      model.N                    = model.RB_stop_Nmax;
      model.ei_numintervals      = [2,2];
      model.RB_numintervals      = [2,2];
    else
      error(['selected model type "', params.model_type, '" is unknown.']);
    end
  end
  model.model_type = params.model_type;

  model = model_default(model);

  model.verbose = 1;
  model.debug   = 0;

end
% vim: set et sw=2:
%| \docupdate 
