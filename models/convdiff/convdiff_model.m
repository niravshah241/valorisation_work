function model=convdiff_model(dummy)
% function model=convdiff_model(dummy)
% function creating a simple model for a linear convection diffusion problem
%
% This model describes a linear fuelcell problem with three parameters.

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


  model = lin_evol_model_default;

  model.name = 'convdiff_model';
  model.rb_problem_type = 'lin_evol';
  %% data that is definitely used outside of the detailed simulations
  model.mu_names           = {'c_init','beta','k'};
  model.mu_ranges          = {[0 1],[0 1],[0 5e-8]};

  %% data that might be of interest outside of the detailed simulations
  model.T             = 0.5;
  model.nt            = 200;
  %model.dt;
  %model.mu;

  %model.N;
  %model.Nmax;
  %model.RB_numintervals;

  model.xrange        = [0,1000e-6]; % grid x-coordinate range
  model.yrange        = [0,200e-6];  % grid y-coordinate range
  model.xnumintervals = 200;  % number of grid cellls in x-direction
  model.ynumintervals = 40;   % number of grid cellls in y-direction

  model.operators_ptr              = @fv_operators_implicit_explicit;
  model.operators_diff_implicit    = @fv_operators_diff_implicit_gradient;
  model.operators_conv_explicit    = @fv_operators_conv_explicit_lax_friedrichs;
  model.operators_neumann_explicit = @fv_operators_neumann_explicit;

  % init values algorithm
  model.init_values_algorithm      = @fv_init_values;
  %  model.init_values_algorithm = @(model,detailed_data) ...
  %    model.x0_function_ptr(model,detailed_data);

  %% specify initial data function
  %    name_init_values = 'homogeneous';
  %    c_init = 0.0;
  model.init_values_ptr = @init_values_wave;
  model.freq_x = 5*2*pi/1000e-6;
  model.freq_y = 0;
  model.c_init = 1.0;
  %    name_init_values = 'homogeneous';
  %    c_init = 0.3;

  %    name_init_values = 'leftright';
  %    c_init_left = 0.3;
  %    c_init_right = 0.7;
  %    init_middle = (params.xrange(2)+params.xrange(1))/2;

  %% convective flux specification
  %    conv_flux_ptr = @conv_flux_gdl2;
  model.conv_flux_ptr = @conv_flux_velocity_matrixfile;
  model.lxf_lambda              = 1.0194e+003;
  model.divclean_mode           = false;
  model.flux_quad_degree        = 1;
  model.use_velocity_matrixfile = 1;

  %model.velocity_matrixfile;
  model.filecache_velocity_matrixfile_extract = 0;

  %% data that is relevant for detailed simulations only

  %% specify grid type
  % generate cartesian grid
  model.gridtype = 'rectgrid';

  %% specify dirichlet boundary values
  model.dirichlet_values_ptr = @dirichlet_values_weighted_boxes;
  model.c_dir                 = 1.0; % maximum value
  model.beta                  = 0.0;  % weighting parameter: 1=box1, 0 = box2

  %% specify neuman boundary values
  model.neumann_values_ptr = @neumann_values_rightflow;
  %    name_neuman_values = 'zero';
  %    name_neuman_values = 'homogeneous';
  %    c_neu = 0;

  %% specify further data functions
  model.flux_linear      = 1;
  model.diffusivity_ptr  = @diffusivity_homogeneous;
  %    k                       = 0.000000;   % diffusion coefficient
  %    k                       = 0.000001;   % diffusion coefficient
  %    k                       = 1e-15;      % diffusion coefficient
  %    k                       = 1e-10;      % diffusion coefficient
  model.k                     = 1e-8;       % diffusion coefficient

  % precomputed, divcleaned velocity field
  model.name_flux               = 'gdl2';

  model.lambda = 2.0e-7;   % v = - lambda * grad p
  model.inner_product_matrix_algorithm = @fv_inner_product_matrix;
  model.mass_matrix          = @fv_mass_matrix;
  model.energy_norm_gamma    = 0.5;
  model.alpha_over_k         = 1.9915e7;
  model.coercivity_bound_ptr = @fv_coercivity_bound;


  % fields to ignore, such that filecaching works in basis generation
  model.filecache_ignore_fields_in_model = {'N','Nmax'};
  model.filecache_ignore_fields_in_detailed_data = {'RB_info'};

  model.RB_generation_mode = 'greedy_uniform_fixed';
  % choose extension method
  model.RB_extension_algorithm = @RB_extension_PCA_fixspace;
  %                 'RB_extension_max_error_snapshot'};
  model.RB_stop_timeout = 60*60; % 1 hour
  model.RB_stop_Nmax = 20;
  model.RB_stop_max_val_train_ratio = inf;

  % target accuracy epsilon:

  model.RB_stop_epsilon = 1e-9; % 0.003 is already realized by init data!!
  %params.stop_on_val_train_ratio = 0;

  % decide, whether estimator or true error is error-indicator for greedy
  % params.RB_error_indicator = 'error'; % true error
  model.RB_error_indicator = 'estimator'; % Delta from rb_simulation

  % test parameters
  model.RB_detailed_test_savepath = 'test_data_100';
  model.RB_test_size = 1000;

  %bnd_rect_corner1;
  %bnd_rect_corner2;
  %bnd_rect_index   = [];

  model.bnd_rect_corner1 = [ model.xrange(1) ,         model.yrange(2) ; ...
                             model.xrange(1)+ 250e-6 , model.yrange(2); ...
                             model.xrange(1) ,         model.yrange(1); ...
                             model.xrange(1) ,         model.yrange(1);
                             model.xrange(2) ,         model.yrange(1)
                           ];
  model.bnd_rect_corner1 = model.bnd_rect_corner1' - eps;
  model.bnd_rect_corner2 = [ model.xrange(2),        model.yrange(2); ...
                             model.xrange(2)-250e-6, model.yrange(2); ...
                             model.xrange(2),        model.yrange(1);
                             model.xrange(1),        model.yrange(2);
                             model.xrange(2),        model.yrange(2)
                           ];
  model.bnd_rect_corner2 = model.bnd_rect_corner2' + eps;
  model.bnd_rect_index   = [-1, -2, -2 , -1, -2 ...
                           ];

  model.use_velocity_matrixfile = 1;

  %in case of later matrix generation: downscale original velocities
  %    model.divclean_downscale = 1;
  %    model.divclean_downscale_quota = 0.04;

  % set matrix-file name for optional generation or reading of file
  model.velocity_matrixfile = [ 'vel_', model.name_flux,'_',...
                                num2str( model.xnumintervals ),'x',...
                                num2str( model.ynumintervals ),...
                                '_l',num2str(    model.lambda),'.mat'];

  model.dir_box_xrange = {[ model.xrange(1)-eps, ...
                            0.5*( model.xrange(1) + model.xrange(2) )], ...
                          [ 0.5*( model.xrange(1) + model.xrange(2) )-eps, ...
                            model.xrange(2)]...
                         };
  model.dir_box_yrange = {[ model.yrange(2)-eps,     model.yrange(2)+eps], ...
                          [ model.yrange(2)-eps,     model.yrange(2)+eps] ...
                         };

  model.RB_numintervals = 4 * ones(size(model.mu_names));

  model = model_default(model);

  model.save_time_indices = 0:model.nt;

  model.verbose = 5;
  model.debug   = 0;
