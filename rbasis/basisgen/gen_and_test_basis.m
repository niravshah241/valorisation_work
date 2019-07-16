function gen_and_test_basis(compute_RB_methods,test_RB_methods)
%function gen_and_test_basis(compute_RB_methods,test_RB_methods)
%
% function performing basis generation by different methods or
% testing of bases. The cell-array compute_RB_methods contains a 
% list of basis-names to generate and save. Similarly the
% test_RB_methods contains the names of the files to read. These
% are extended and stored in a new filename extended by '_tested'.

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


% Bernard Haasdonk 9.6.2007

%%% basic structures

% make sure, that the following path exists!
basisgen_savepath = fullfile(getenv('RBMATLABTEMP'),'basisgen');
if ~exist(basisgen_savepath)
  error('please create savepath!');
end;

% select here, which methods of the ones below are to be called 
% for construction

%compute_RB_methods = {}

% for paper-experiments:
%compute_RB_methods = {...
%    'RB_adaptive_refined_2_2_r1_smax1000_theta5e-2',...
%    'RB_adaptive_refined_2_2_r5_smax1000_theta5e-2',...
%    'RB_adaptive_refined_2_2_r1_smax1000_theta1e-1',...
%    'RB_adaptive_refined_2_2_r5_smax1000_theta1e-1',...
%    'RB_uniform_refined_2_2_r5',...
%		   };
%compute_RB_methods = {...
%    'RB_uniform_refined_2_2_r1',...
%    'RB_uniform_refined_2_2_r5',...
%    'RB_uniform_fixed_64',...
%    'RB_uniform_fixed_256',...
%    'RB_uniform_refined_2_2_r1',...
%    'RB_uniform_refined_2_2_r5',...
%    'RB_adaptive_refined_2_2_r1_smax1000_theta5e-2',...
%    'RB_adaptive_refined_2_2_r5_smax1000_theta5e-2',...
%    'RB_adaptive_refined_2_2_r1_smax1000_theta1e-1',...
%    'RB_adaptive_refined_2_2_r5_smax1000_theta1e-1',...
%		   };




%compute_RB_methods = {...
%    'RB_adaptive_refined_2_2_r1_smax1000',...
%    'RB_adaptive_refined_2_2_r1_smax100'}; %,...

%    'RB_adaptive_refined_2_2_r1_smax10',...
%    'RB_adaptive_refined_2_2_r1_smax5',...
%    'RB_adaptive_refined_2_2_r1_smax1'};

%compute_RB_methods = { 'RB_adaptive_refined_3_3',...
%		    'RB_uniform_refined_3_3',...
%		    'RB_random_fixed_50',...
%		    'RB_random_fixed_100',...
%		    'RB_uniform_fixed_49',...
%		    'RB_uniform_fixed_100' ...
%		   };
%compute_RB_methods = {'RB_random_fixed_100'};
%compute_RB_methods = {'RB_uniform_fixed_49'};
%compute_RB_methods = {'RB_uniform_refined_3_3'};
%compute_RB_methods = {'RB_adaptive_refined_3_3'};
%compute_RB_methods = {'RB_adaptive_refined_3_3_Nmax40'};
%compute_RB_methods = {'RB_adaptive_refined_3_3'};
%compute_RB_methods = {'RB_adaptive_refined_2_2',...
%		      'RB_uniform_fixed_100' ...
%		   };

%compute_RB_methods = {'RB_adaptive_refined_2_2_r1', ...
%		    'RB_adaptive_refined_2_2_r2', ...
%		    'RB_adaptive_refined_2_2_r5', ...
%		    'RB_adaptive_refined_2_2_r10'};
%compute_RB_methods = {'RB_uniform_fixed_16', ...
%		    'RB_uniform_fixed_64', ...
%		    'RB_uniform_fixed_256', ...
%		    'RB_uniform_fixed_1024', ...
%		    'RB_random_fixed_16', ...
%		    'RB_random_fixed_64', ...
%		    'RB_random_fixed_256', ...
%		    'RB_random_fixed_1024'};

%compute_RB_methods = {};

% select here, which methods of the ones below are to be called 
% for test

%test_RB_methods = compute_RB_methods;

if nargin<2 
  test_RB_methods = {};
end;

%test_RB_methods = {...
%    'RB_uniform_refined_2_2_r1',...
%    'RB_uniform_refined_2_2_r5',...
%    'RB_uniform_fixed_64',...
%    'RB_uniform_fixed_256',...
%    'RB_uniform_refined_2_2_r1',...
%    'RB_uniform_refined_2_2_r5',...
%    'RB_adaptive_refined_2_2_r1_smax1000_theta5e-2',...
%    'RB_adaptive_refined_2_2_r5_smax1000_theta5e-2',...
%    'RB_adaptive_refined_2_2_r1_smax1000_theta1e-1',...
%    'RB_adaptive_refined_2_2_r5_smax1000_theta1e-1',...
%		   };

%test_RB_indicators = {}; % none
%test_RB_indicators = {'error','estimator'};
test_RB_indicators = {'estimator'}; % only estimators

% set initial parameters
%orgparams = basisgen_init_params;

disp('1. setting common parameters.');

params = [];

params.xrange = [0,1000e-6]; % grid x-coordinate range
params.yrange = [0,200e-6];  % grid y-coordinate range
params.T = 0.5;              % maximum time
params.xnumintervals = 200;             % number of grid cellls in x-direction
params.ynumintervals = 40;              % number of grid cellls in y-direction
params.nt = 200;             % number of time steps 

% set upper boundary to dirichlet, middle of upper and the lower boundary
% to neuman
% left set to dirichlet, right set to outflow

params.bnd_rect_corner1 = [params.xrange(1) , params.yrange(2) ; ...
		    params.xrange(1)+ 250e-6 , params.yrange(2); ...
		    params.xrange(1) , params.yrange(1); ...
		    params.xrange(1) , params.yrange(1);
		    params.xrange(2) , params.yrange(1)
		   ];
params.bnd_rect_corner1 = params.bnd_rect_corner1' - eps;
params.bnd_rect_corner2 = [params.xrange(2), params.yrange(2); ...
		    params.xrange(2)-250e-6, params.yrange(2); ...
		    params.xrange(2), params.yrange(1);
		    params.xrange(1),params.yrange(2);
		    params.xrange(2),params.yrange(2)
		   ];
params.bnd_rect_corner2 = params.bnd_rect_corner2' + eps;

params.bnd_rect_index = [-1, -2, -2 , -1, -2 ...
		   ];
%params.name_dirichlet_values = 'homogeneous';
%params.c_dir = 1.0;

% nur oben nonzero values, darunter 0 von links hinein.
%params.name_dirichlet_values = 'uplow';
%params.c_dir_up = 1.0;
%params.c_dir_low = 000.0;
%params.dir_middle = params.yrange(2)-eps;

% nur Einfluss linke obere Ecke
%params.name_dirichlet_values = 'box';
%params.dir_box_xrange = [params.xrange(1)-eps, 0.5*(params.xrange(1) ...
%						  + params.xrange(2))];
%params.dir_box_yrange = [params.yrange(2)-eps, params.yrange(2)+eps];
%params.c_dir = 1.0;

% nur Einfluss linke obere Ecke oder rechte obere (Konvexcombination)
params.name_dirichlet_values = 'weighted_boxes';
params.dir_box_xrange = {[params.xrange(1)-eps, ...
		    0.5*(params.xrange(1)+params.xrange(2))], ...
		    [0.5*(params.xrange(1)+params.xrange(2))-eps, ...
		    params.xrange(2)]...
		   };
params.dir_box_yrange = {[params.yrange(2)-eps, params.yrange(2)+eps], ...
		    [params.yrange(2)-eps, params.yrange(2)+eps] ...
		   };
params.c_dir = 1.0; % maximum value
params.beta = 0.0; % weighting parameter: 1=box1, 0 = box2 

params.name_neuman_values = 'rightflow';
%params.name_neuman_values = 'zero';
%params.name_neuman_values = 'homogeneous';
%params.c_neu = 0;

% parameters inducing "affine dependency" into the problem
params.name_diffusivity = 'homogeneous';
%params.name_diffusive_num_flux = 'gradient2';
params.name_diffusive_num_flux = 'gradient';
%params.k = 0.000000;      % diffusion coefficient 
%params.k = 0.000001;      % diffusion coefficient 
%params.k = 1e-15;      % diffusion coefficient 
%params.k = 1e-10;      % diffusion coefficient 
params.k = 1e-5;      % diffusion coefficient 
		  
%params.name_flux = 'gdl_pgrad_and_parabola';		       
%params.v_in_max = 0.0;          % maximum inflow velocity 0 to 1 is OK

% precomputed, divcleaned velocity field
params.name_flux = 'gdl2';		       
% if precomputed div-cleaned flux is wanted, set this flag.
disp('using precomputed velocity matrix file!!');
params.use_velocity_matrixfile = 1;
params.divclean_mode = 'none'; % file is already cleaned by optimization

%in case of later matrix generation: downscale original velocities
%params.divclean_downscale = 1;
%params.divclean_downscale_quota = 0.04;

params.lambda = 2.0e-7;   % v = - lambda * grad p
% set matrix-file name for optional generation or reading of file 
params.velocity_matrixfile = ['vel_',params.name_flux,'_',...
		    num2str(params.xnumintervals),'x',...
		    num2str(params.ynumintervals),...
		    '_l',num2str(params.lambda),'.mat'];

params.name_convective_num_flux = 'lax-friedrichs';
%params.lxf_lambda = 0.1172; % lambda value for Lax-Friedrichs diffusivity
%params.lxf_lambda = 0.001172; % lambda value for Lax-Friedrichs diffusivity

%params.name_init_values = 'homogeneous';
%params.c_init = 0.0;

params.name_init_values = 'wave';
params.freq_x = 5*2*pi/1000e-6;
params.freq_y = 0;
params.c_init = 1.0;

%params.name_init_values = 'homogeneous';
%params.c_init = 0.3;

%params.name_init_values = 'leftright';
%params.c_init_left = 0.3;
%params.c_init_right = 0.7;
%params.init_middle = (params.xrange(2)+params.xrange(1))/2;

% output options:
params.verbose = 10;  % 0   : no output, 
                      % >9  : informative output
                      % >19 : verbose informative output
                      % >29 : debug output 
                      % >99 : desperate debug output, with dbstops
params.show_colorbar = 1;
params.colorbar_mode = 'EastOutside';                 

% generate cartesian grid
params.gridtype =  'rectgrid'; 

%params.epsilon = 1e-6; 
%params.epsilon = 1e-4; % => No refinement required, single extension loop..!!
%params.affine_decomp_mode= 'none';
params.L_I_inv_norm_bound = 1;
params.L_E_norm_bound = 1;
%params.detailed_simulation_algorithm = 'fv_conv_diff_operators';
params.detailed_simulation_algorithm = 'detailed_simulation';
params.operators_algorithm = 'fv_operators_implicit_explicit';
params.init_values_algorithm = 'fv_init_values';
params.inner_product_matrix_algorithm = 'fv_inner_product_matrix';
params.error_algorithm = 'fv_error';
params.lxf_lambda = 1.0194e+003;
params.data_const_in_time = 1;
params.error_norm = 'l2';
%keyboard;

% set rb-parameter information: 2d, cinit fixed to 1
% in 3d-parameter set examples below ("p3" stringpart) these
% settings are changed.
params.mu_names = {'beta','k'};
params.mu_ranges = {[0 1],[0 5e-8]};

% choose extension method
params.RB_extension_algorithm = 'RB_extension_PCA_fixspace';
% 		    'RB_extension_max_error_snapshot'}; 
params.RB_stop_timeout = 60*60; % 1 hour 				
				
% target accuracy epsilon:

params.RB_stop_epsilon = 1e-9; % 0.003 is already realized by init data!!
%params.stop_on_val_train_ratio = 0;

% decide, whether estimator or true error is error-indicator for greedy
% params.RB_error_indicator = 'error'; % true error
params.RB_error_indicator = 'estimator'; % Delta from rb_simulation

% test parameters
params.RB_detailed_test_savepath = 'test_data_100';
params.RB_test_size = 1000;

params.rb_problem_type = 'lin_evol';

% --------------------------------------------------------------------------

% store original params

orgparams = params;

for meth = 1:length(compute_RB_methods)

  disp('--------------------------------------------------------------')
  disp(['starting basis construction for ',compute_RB_methods{meth}]);
  disp('clearing cache for fair time comparisons')
  filecache_clear;
  % reload params, such that independent of order of loop
  params = orgparams;

  switch compute_RB_methods{meth}
    
    % ------------------------------------- --------------------------
    % ------------------- 2D Parameterspace --------------------------
    % ------------------------------------- --------------------------
    
   case {'RB_random_fixed_16', ...
	 'RB_random_fixed_64', ...
	 'RB_random_fixed_256', ...
	 'RB_random_fixed_1024'}
    
    params.RB_stop_Nmax = 130; 
    rpos = findstr('_',compute_RB_methods{meth});
    ntrain = str2num(compute_RB_methods{meth}(rpos(end)+1:end));
    disp(['extracted ntrain value = ',num2str(ntrain)]);
    
    params.RB_generation_mode = 'random_fixed'
    params.RB_train_size = ntrain;
    params.RB_train_rand_seed = 12345;
    detailed_data = rb_detailed_prep(params);

   case {'RB_uniform_fixed_16', ...
	 'RB_uniform_fixed_64', ...
	 'RB_uniform_fixed_256', ...
	 'RB_uniform_fixed_1024'}
    
    params.RB_stop_Nmax = 130; 
    rpos = findstr('_',compute_RB_methods{meth});
    ntrain = str2num(compute_RB_methods{meth}(rpos(end)+1:end));
    numintervals = sqrt(ntrain)-1;
    disp(['extracted numintervals value = ',num2str(numintervals)]);
    
    %params.RB_extension_stop = 'epsilon';
    params.RB_generation_mode = 'uniform_fixed'
    params.RB_numintervals = ones(1,2) * numintervals;
    detailed_data = rb_detailed_prep(params);    
    
   case 'RB_adaptive_refined_3_3'%% 3. adaptive refined method:
    
    params.RB_stop_Nmax = 130; 
    params.RB_generation_mode = 'refined';
    params.RB_refinement_mode = 'adaptive';
    params.RB_element_indicator_mode = 'nodes_cogs';
    params.RB_numintervals = [2,2];
    %params.RB_extension_stop = 'epsilon_or_pe_ratio';
    %params.epsilon = 5e-5;
    %params.epsilon = 1e-3;
    %    params.RB_stop_max_val_train_ratio = 1.2;
    params.RB_stop_max_val_train_ratio = 1.2;
    params.RB_M_val_size = 10;
    params.RB_val_rand_seed = 2345;
    params.RB_refinement_theta = 0.05;
    params.RB_max_refinement_level = 15;
    
    detailed_data = rb_detailed_prep(params);
    
   case 'RB_adaptive_refined_2_2'%% 3. adaptive refined method:
    
    params.RB_stop_Nmax = 130; 
    params.RB_generation_mode = 'refined';
    params.RB_refinement_mode = 'adaptive';
    params.RB_element_indicator_mode = 'nodes_cogs';
    params.RB_numintervals = [1,1];
    %params.RB_extension_stop = 'epsilon_or_pe_ratio';
    %params.epsilon = 5e-5;
    %params.epsilon = 1e-3;
    %    params.RB_stop_max_val_train_ratio = 1.2;
    params.RB_stop_max_val_train_ratio = 1.2;
    params.RB_M_val_size = 10;
    params.RB_val_rand_seed = 2345;
    params.RB_refinement_theta = 0.05;
    params.RB_max_refinement_level = 8;
    
    detailed_data = rb_detailed_prep(params);
    
   case {'RB_adaptive_refined_2_2_r1', ...
	 'RB_adaptive_refined_2_2_r2', ...
	 'RB_adaptive_refined_2_2_r5', ...
	 'RB_adaptive_refined_2_2_r10'} %% 3. adaptive refined method:
    
    params.RB_stop_Nmax = 130; 
    rpos = findstr('r',compute_RB_methods{meth});
    rvalue = str2num(compute_RB_methods{meth}(rpos(end)+1:end));
    disp(['extracted r value = ',num2str(rvalue)]);
    
    params.RB_generation_mode = 'refined';
    params.RB_refinement_mode = 'adaptive';
    params.RB_element_indicator_mode = 'nodes_skippedrefs';
    params.RB_element_indicator_s_max = 10;
    
    params.RB_numintervals = [1,1];
    %params.RB_extension_stop = 'epsilon_or_pe_ratio';
    %params.epsilon = 5e-5;
    %params.epsilon = 1e-3;
    %    params.RB_stop_max_val_train_ratio = 1.2;
    params.RB_stop_max_val_train_ratio = rvalue;
    params.RB_M_val_size = 10;
    params.RB_val_rand_seed = 2345;
    params.RB_refinement_theta = 0.05;
    params.RB_max_refinement_level = 9;

    detailed_data = rb_detailed_prep(params);
   case {'RB_adaptive_refined_2_2_r1_smax10',...
	 'RB_adaptive_refined_2_2_r1_smax1',...
	 'RB_adaptive_refined_2_2_r1_smax5',...
	 'RB_adaptive_refined_2_2_r1_smax1000',...
	 'RB_adaptive_refined_2_2_r1_smax100'}; 
    
    params.RB_stop_Nmax = 130; 
    rpos = findstr('x',compute_RB_methods{meth});
    smax = str2num(compute_RB_methods{meth}(rpos(end)+1:end));
    disp(['extracted smax value = ',num2str(smax)]);
    
    params.RB_generation_mode = 'refined';
    params.RB_refinement_mode = 'adaptive';
    params.RB_element_indicator_mode = 'nodes_skippedrefs';
    params.RB_element_indicator_s_max = smax;
    
    params.RB_numintervals = [1,1];
    %params.RB_extension_stop = 'epsilon_or_pe_ratio';
    %params.epsilon = 5e-5;
    %params.epsilon = 1e-3;
    %    params.RB_stop_max_val_train_ratio = 1.2;
    params.RB_stop_max_val_train_ratio = 1;
    params.RB_M_val_size = 10;
    params.RB_val_rand_seed = 2345;
    params.RB_refinement_theta = 0.05;
    params.RB_max_refinement_level = 9;

    detailed_data = rb_detailed_prep(params);

   case { 'RB_adaptive_refined_2_2_r1_smax1000_theta5e-2',...
	  'RB_adaptive_refined_2_2_r5_smax1000_theta5e-2',...
	  'RB_adaptive_refined_2_2_r1_smax1000_theta1e-1',...
	  'RB_adaptive_refined_2_2_r5_smax1000_theta1e-1'}
    
    params.RB_stop_Nmax = 130; 
    smaxpos = findstr('smax',compute_RB_methods{meth});
    rpos = findstr('r',compute_RB_methods{meth});
    thetapos = findstr('theta',compute_RB_methods{meth});
    rvalue = str2num(compute_RB_methods{meth}(rpos(end)+1:smaxpos(end)-2));
    disp(['extracted rvalue = ',num2str(rvalue)]);
    smax = str2num(compute_RB_methods{meth}(smaxpos(end)+4:thetapos(end)-2));
    disp(['extracted smax value = ',num2str(smax)]);
    theta = str2num(compute_RB_methods{meth}(thetapos(end)+5:end));
    disp(['extracted theta value = ',num2str(theta)]);
    
    params.RB_generation_mode = 'refined';
    params.RB_refinement_mode = 'adaptive';
    params.RB_element_indicator_mode = 'nodes_skippedrefs';
    params.RB_element_indicator_s_max = smax;
    
    params.RB_numintervals = [1,1];
    %params.RB_extension_stop = 'epsilon_or_pe_ratio';
    %params.epsilon = 5e-5;
    %params.epsilon = 1e-3;
    %    params.RB_stop_max_val_train_ratio = 1.2;
    params.RB_stop_max_val_train_ratio = rvalue;
    params.RB_M_val_size = 10;
    params.RB_val_rand_seed = 2345;
    params.RB_refinement_theta = theta;
    params.RB_max_refinement_level = 9;

    detailed_data = rb_detailed_prep(params);

    
   case 'RB_adaptive_refined_3_3_Nmax40'

    params.RB_stop_Nmax = 130; 
    params.RB_generation_mode = 'refined';
    params.RB_refinement_mode = 'adaptive';
    params.RB_element_indicator_mode = 'nodes_cogs';
    params.RB_numintervals = [2,2];
    %params.RB_extension_stop = 'epsilon_or_pe_ratio';
    %params.epsilon = 5e-5;
    %params.epsilon = 1e-3;
    %    params.RB_stop_max_val_train_ratio = 1.2;
    params.RB_stop_max_val_train_ratio = 1.2;
    params.RB_M_val_size = 10;
    params.RB_val_rand_seed = 2345;
    params.RB_refinement_theta = 0.05;
    params.RB_max_refinement_level = 15;
    params.RB_stop_Nmax = 40;
    
    detailed_data = rb_detailed_prep(params);
   
   case {'RB_uniform_refined_2_2_r1',...   
	 'RB_uniform_refined_2_2_r5'}
    
    params.RB_stop_Nmax = 130; 
    rpos = findstr('r',compute_RB_methods{meth});
    rvalue = str2num(compute_RB_methods{meth}(rpos(end)+1:end));
    disp(['extracted r value = ',num2str(rvalue)]);
    
    params.RB_generation_mode = 'refined';
    params.RB_refinement_mode = 'uniform';
    params.RB_numintervals = [2,2];
    %    params.RB_stop_max_val_train_ratio = 1.2;
    params.RB_stop_max_val_train_ratio = rvalue;
    params.RB_M_val_size = 10;
    params.RB_val_rand_seed = 2345;
    params.RB_max_refinement_level = 15;
    
    detailed_data = rb_detailed_prep(params);

   case 'RB_uniform_refined_3_3'   
    %% 4. uniform refined method
    
    params.RB_stop_Nmax = 130; 
    params.RB_generation_mode = 'refined';
    params.RB_refinement_mode = 'uniform';
    params.RB_numintervals = [2,2];
%    params.RB_stop_max_val_train_ratio = 1.2;
    params.RB_stop_max_val_train_ratio = 1.2;
    params.RB_M_val_size = 10;
    params.RB_val_rand_seed = 2345;
    params.RB_max_refinement_level = 15;
    
    detailed_data = rb_detailed_prep(params);

    % ------------------------------------- --------------------------
    % ------------------- 3D Parameterspace --------------------------
    % ------------------------------------- --------------------------
    
    
   case {'RB_uniform_fixed_p3_27', ...
	 'RB_uniform_fixed_p3_64', ...
	 'RB_uniform_fixed_p3_125'}
    

    % set rb-parameter information: 3d
    params.mu_names = {'c_init','beta','k'};
    params.mu_ranges = {[0,1],[0 1],[0 5e-8]};
    params.RB_stop_Nmax = 150; 
    
    rpos = findstr('_',compute_RB_methods{meth});
    ntrain = str2num(compute_RB_methods{meth}(rpos(end)+1:end));
    numintervals = round(ntrain^(1/3))-1;
    disp(['extracted numintervals value = ',num2str(numintervals)]);
    
    %params.RB_extension_stop = 'epsilon';
    params.RB_generation_mode = 'uniform_fixed'
    params.RB_numintervals = ones(1,3) * numintervals;
    detailed_data = rb_detailed_prep(params);    

   case {'RB_uniform_fixed_p3_216_notimeout'}

    % set rb-parameter information: 3d
    params.mu_names = {'c_init','beta','k'};
    params.mu_ranges = {[0,1],[0 1],[0 5e-8]};
    params.RB_stop_Nmax = 150; 
    params.RB_stop_timeout = inf;
    
    rpos = findstr('_',compute_RB_methods{meth});
    tpos = findstr('notimeout',compute_RB_methods{meth});
    
    ntrain = str2num(compute_RB_methods{meth}(rpos(end-1)+1:tpos(end)-2));
    numintervals = round(ntrain^(1/3))-1;
    disp(['extracted numintervals value = ',num2str(numintervals)]);
    
    %params.RB_extension_stop = 'epsilon';
    params.RB_generation_mode = 'uniform_fixed'
    params.RB_numintervals = ones(1,3) * numintervals;
    detailed_data = rb_detailed_prep(params);    
    
   case { ...
       'RB_adaptive_refined_p3_i2_r1_smax1000_theta5e-2',...
       'RB_adaptive_refined_p3_i3_r1_smax1000_theta5e-2',...
       'RB_adaptive_refined_p3_i3_r1_smax1000_theta5e-2_notimeout'...
	};
    
    params.mu_names = {'c_init','beta','k'};
    params.mu_ranges = {[0,1],[0 1],[0 5e-8]};
    params.RB_stop_Nmax = 150; 
        
    ipos = findstr('_i',compute_RB_methods{meth});
    smaxpos = findstr('smax',compute_RB_methods{meth});
    rpos = findstr('r',compute_RB_methods{meth});
    thetapos = findstr('theta',compute_RB_methods{meth});
    tpos = findstr('notimeout',compute_RB_methods{meth});
    ivalue = str2num(compute_RB_methods{meth}(ipos(end)+2:rpos(end)-2));
    disp(['extracted ivalue = ',num2str(ivalue)]);
    rvalue = str2num(compute_RB_methods{meth}(rpos(end)+1:smaxpos(end)-2));
    disp(['extracted rvalue = ',num2str(rvalue)]);
    smax = str2num(compute_RB_methods{meth}(smaxpos(end)+4:thetapos(end)-2));
    disp(['extracted smax value = ',num2str(smax)]);
    if isempty(tpos)
      theta = str2num(compute_RB_methods{meth}(thetapos(end)+5:end));
    else      
      theta = str2num(compute_RB_methods{meth}(thetapos(end)+5:tpos(end)-2));
      params.RB_stop_timeout = inf;
    end;
    disp(['extracted theta value = ',num2str(theta)]);
    
    params.RB_generation_mode = 'refined';
    params.RB_refinement_mode = 'adaptive';
    params.RB_element_indicator_mode = 'nodes_cogs_skippedrefs';
    params.RB_element_indicator_s_max = smax;
    
    params.RB_numintervals = ones(1,3)*ivalue;
    %params.RB_extension_stop = 'epsilon_or_pe_ratio';
    %params.epsilon = 5e-5;
    %params.epsilon = 1e-3;
    %    params.RB_stop_max_val_train_ratio = 1.2;
    params.RB_stop_max_val_train_ratio = rvalue;
    params.RB_M_val_size = 10;
    params.RB_val_rand_seed = 2345;
    params.RB_refinement_theta = theta;
    params.RB_max_refinement_level = 15;
    
    detailed_data = rb_detailed_prep(params);
    
   case {...
       'RB_uniform_refined_p3_i2_r1',...
       'RB_uniform_refined_p3_i3_r1',...
       'RB_uniform_refined_p3_i2_r1_notimeout',...
       'RB_uniform_refined_p3_i3_r1_notimeout'...
	}
    
    params.mu_names = {'c_init','beta','k'};
    params.mu_ranges = {[0,1],[0 1],[0 5e-8]};
    params.RB_stop_Nmax = 150; 
    
    ipos = findstr('_i',compute_RB_methods{meth});
    rpos = findstr('r',compute_RB_methods{meth});
    tpos = findstr('notimeout',compute_RB_methods{meth});
    ivalue = str2num(compute_RB_methods{meth}(ipos(end)+2:rpos(end)-2));
    disp(['extracted ivalue = ',num2str(ivalue)]);
    if isempty(tpos)
      rvalue = str2num(compute_RB_methods{meth}(rpos(end)+1:end));
    else
      rvalue = str2num(compute_RB_methods{meth}(rpos(end)+1:tpos(end)-2));
      params.RB_stop_timeout = inf;
    end;
    
    disp(['extracted r value = ',num2str(rvalue)]);
    
    params.RB_generation_mode = 'refined';
    params.RB_refinement_mode = 'uniform';
    params.RB_numintervals = ones(1,3)*ivalue;
    %    params.RB_stop_max_val_train_ratio = 1.2;
    params.RB_stop_max_val_train_ratio = rvalue;
    params.RB_M_val_size = 10;
    params.RB_val_rand_seed = 2345;
    params.RB_max_refinement_level = 15;
    
    detailed_data = rb_detailed_prep(params);
    
   otherwiseh
    error('basis-extension method unknown!');
    
  end;
    
  save(fullfile(basisgen_savepath,compute_RB_methods{meth}),...
       'detailed_data','params');
  
end; % loop RB-generation methods

%%% Test: Assuming test data to be out of the synchronize directory
%called data/test_data

if ~isempty(test_RB_methods)
  % generate random parameter vectors data
    
  for meth = 1:length(test_RB_methods)
    
    fn = fullfile(basisgen_savepath,test_RB_methods{meth}); 
    load(fn); 
    
    %    params = orgparams;
    grid =  construct_grid(params); 
    
    rand('state',1111);
    M = rand_uniform(params.RB_test_size, params.mu_ranges);
    
    if ismember(test_RB_indicators,'error');
      save_detailed_simulations(M, params.RB_detailed_test_savepath, ...
				grid, params);  
    end;
    % extract final slices for faster test
    %save_detailed_simulation_finals(params.test_path);  
    
    % assuming, that at least stat and RB is contained!!
    %  [stat_test] = rb_test_basis(RB, stat,params);
    %  params.test_path = fullfile(rbmatlabhome,'tmpdata','test_data_10');
    
    tmpparams = params;
    tmpparams.test_N_samples = 10;
    detailed_data.RB_info.M_test = M;
    
    if ismember(test_RB_indicators,'error');      
      tmpparams.RB_error_indicator = 'error';
      disp(['testing ', fn,' by assessing true errors']);
      [maxerr, max_mu, minerr, min_mu] = ...
	  rb_test_convergence(detailed_data,tmpparams);
      
      detailed_data.RB_info.max_test_error_sequence = maxerr;
      detailed_data.RB_info.min_test_error_sequence = minerr;
    end;
    
    if ismember(test_RB_indicators,'estimator');      
      tmpparams.RB_error_indicator = 'estimator';
      disp(['testing ', fn,' by assessing posterior-estimators']);
      [maxerr, max_mu, minerr,min_mu] = ...
	  rb_test_convergence(detailed_data,tmpparams);
      
      detailed_data.RB_info.max_test_estimator_sequence = maxerr;
      detailed_data.RB_info.min_test_estimator_sequence = minerr;
    end;
    
    save([fn,'_tested'],'detailed_data','params');
  end;
  
end; % end test loop


% TO BE ADJUSTED TO NEW SYNTAX
%| \docupdate 
