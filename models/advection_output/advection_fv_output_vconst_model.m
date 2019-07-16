function model = advection_fv_output_vconst_model(params);
%function model = advection_fv_output_vconst_model(params);
%
% model for the new advection model with functional output
% two overlapping velocity fields and boundary value modification
% are the parameters. Output functional is average concentration 
% in subdomain. Pure explicit discretization of convection with fv
% discretization. In contrast to the advection_fv_output_model, 
% the velocity field is constant in time (but parameter
% dependent). The boundary conditions are decaying in time.
%
% possible fields of params:
%      coarse_factor: coarsening factor between 1 and 8 
%                     time-steps, gridsize are scaled down with this. 
%
% To be used for results of article at-Automatisierungstechnik 2010.

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


% B. Haasdonk 26.8.2009

model = [];
if (nargin ==1) & (isfield(params,'coarse_factor'))
  coarse_factor = params.coarse_factor;
else
  coarse_factor = 1;
end;

% specification of the time information
model.T             = 1;
model.nt            = 2048/coarse_factor;
%model.save_time_indices = 0:16/coarse_factor:model.nt;
model.save_time_indices = 0:8/coarse_factor:model.nt;
disp('nt to be adjusted later!');
model.dt = model.T/model.nt;
model.verbose = 9;
model.debug = 0;
%model.debug = 1; %
%if model.debug
%  disp('model with debugging turned on.');
%end;
model.axis_equal = 1;
model.error_estimation = 1;

% grid information
%model.grid_initfile = 'rectangle_triagrid.mat';
model.gridtype = 'triagrid';
model.xnumintervals = 256/coarse_factor;
model.ynumintervals = 128/coarse_factor;
model.xrange = [0,2];
model.yrange = [0,1];
% set completely dirichlet boundary
model.bnd_rect_corner1 = [0,0]-eps;
model.bnd_rect_corner2 = [2,1]+eps;
%model.opts.bnd_rect_index   = [-1];
model.element_quadrature = @triaquadrature;
model.intersection_quadrature = @intervalquadrature;
model.dim_U = model.xnumintervals * model.ynumintervals * 2;

% Main global function pointers expected in every model
model.gen_model_data =      @lin_evol_gen_model_data;
model.gen_detailed_data =   @lin_evol_gen_detailed_data;
model.gen_reduced_data =     @lin_evol_gen_reduced_data;
model.detailed_simulation = @lin_evol_detailed_simulation;
model.rb_simulation =       @lin_evol_rb_simulation;
model.rb_reconstruction =   @lin_evol_rb_reconstruction;
model.plot_sim_data =       @lin_evol_plot_sim_data;

% Dirichlet and Initial data
model.dirichlet_values_ptr = @dirichlet_values_affine_decomposed;
model.dirichlet_values_coefficients_ptr = @my_dirichlet_values_coefficients;
model.dirichlet_values_components_ptr = @my_dirichlet_values_components;
% dummy one dirichlet_values:
%model.dirichlet_values_coefficients_ptr = @(params) 1;
%model.dirichlet_values_components_ptr = @(glob,params) {ones(1, ...
%						  size(glob,2))};
model.rb_init_data_basis = @RB_init_data_basis;
model.init_values_ptr = @init_values_affine_decomposed;
model.init_values_coefficients_ptr = @my_dirichlet_values_coefficients_t0;
model.init_values_components_ptr = @my_dirichlet_values_components;
model.cone_number = 3; % number of components = cones 
model.cone_range = [0,1]; % x-range of cone-support 
model.cone_weight = 1; % leftmost cone with full weight
model.cone_amplitude = 1; % full amplitude
model.velocity_ptr = @velocity_affine_decomposed;
model.velocity_coefficients_ptr = @my_velocity_coefficients;
model.velocity_components_ptr = @my_velocity_components;
%model.vx_weight = 1.0;
model.vx_weight = 0.75;
%model.vy_weight = 1.0;
model.vy_weight = 0.75;
model.conv_flux_ptr = @conv_flux_linear;

%  set parameter setting function

model.set_time = @set_time_default;
model.set_mu = @set_mu_default;
model.get_mu = @get_mu_default;

% perhaps these are redundant later...
model.divclean_mode = 0;
model.flux_quad_degree = 1;
model.flux_linear = 1;

model.init_values_algorithm = @disc_init_values;
model.init_values_qdeg = 0;
model.pdeg = 0; % FV schemes with piecewise constant ansatz functions
model.evaluate_basis = @fv_evaluate_basis;
model.l2project = @l2project;
model.local_mass_matrix = @fv_local_mass_matrix_tria; % triangular grid
model.mass_matrix = @fv_mass_matrix; 
model.ndofs_per_element = 1;
model.plot = @fv_plot; % then also plot_sequence will work

% for use of old datafunctions, use the following:
%model.operators_ptr = @fv_operators_implicit_explicit;
%model.operators_conv_explicit = @fv_operators_conv_explicit;
%model.operators_conv_implicit = @fv_operators_conv_implicit;
%model.operators_diff_explicit = @fv_operators_diff_explicit;
%model.operators_diff_implicit = @fv_operators_diff_implicit;
%model.operators_neumann_explicit = @fv_operators_neumann_explicit_old;
%model.operators_neumann_implicit = @fv_operators_neumann_implicit;

model.operators_ptr = @fv_operators_implicit_explicit;
% for use of new data function access, use now
model.operators_conv_explicit = @fv_operators_conv_explicit_engquist_osher;
model.operators_conv_implicit = @fv_operators_zero;
model.operators_diff_explicit = @fv_operators_zero;
model.operators_diff_implicit = @fv_operators_zero;
model.operators_neumann_implicit = @fv_operators_zero;
model.operators_neumann_explicit = @fv_operators_zero;
model.operators_output = @fv_operators_output;
%model.operators_diff_implicit = @fv_operators_diff_implicit_gradient;
%model.operators_neumann_explicit = @fv_operators_neumann_explicit;

%model.flux_linear = 1;
model.data_const_in_time = 0; % time varying data...
%model.diffusivity_ptr = ...;
%model.neumann_value_ptr = @...
model.affinely_decomposed = 1; % data functions allow affine
                               % parameter decomposition

model.compute_output_functional = 1; % turn on computation of output
model.output_function_ptr = @output_function_box_mean;
%model.sbox_xmin = 1;
model.sbox_xmin = 1;
model.sbox_xmax = 2;
model.sbox_ymin = 0;
%model.sbox_ymax = 0.5;
model.sbox_ymax = 0.5;

%model.disc_element_mean_ptr = @fv_element_mean;

% RB information:
%model.mu_names           = {'cone_weight','vx_weight','vy_weight'};
model.mu_names           = {'cone_amplitude','vx_weight','vy_weight'};
model.mu_ranges          = {[0 1],[0 1],[0,1]};

return; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% auxiliary functions used as pointers above:
%%%%% check later, if they are of general interest to be exported
%%%%% as standalone functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% new function syntax: global coordinates as rows in glob_coord,
% time and parameter variables in params

% dirichlet-data: convex combination of equidistant cones, decaying
% in time
function res = my_dirichlet_values_coefficients_t0(params);
params.t = 0;
res = my_dirichlet_values_coefficients(params);

function res = my_dirichlet_values_coefficients(params);
% res is a vector of coefficients
Q_0 = params.cone_number;
res = zeros(1,Q_0);
max_pos = params.cone_weight * (Q_0-1)+1;
t = params.t;
for q = 1:Q_0  
  res(q) = (1-(max_pos-q))*(1-t) * ((max_pos>=q) && (max_pos < q+1)) ...
	   + (1+(max_pos-q))*(1-t) * ((max_pos>=q-1) && (max_pos < q));
end;
res = res * params.cone_amplitude;

function res = my_dirichlet_values_components(glob,params);
% res is a cell-array of row-vectors of global evaluations
Q_0 = params.cone_number;
res = cell(1,Q_0);
%delta_cone = 1/(Q_0+1);
delta_cone = (params.cone_range(2)-params.cone_range(1))/(Q_0+1);
cone_pos_x = delta_cone * (1:Q_0)+ params.cone_range(1);
for q = 1:Q_0  
  res{q} = 1-min(sqrt((glob(:,1)-cone_pos_x(q)).^2+...
	    (glob(:,2)-1).^2)*(Q_0+1),1);
end;

% velocity field: overlap of y- and x parabolic profiles
function res = my_velocity_coefficients(params);
%res = [params.vx_weight, params.vy_weight]*(1-params.t);
res = [params.vx_weight, params.vy_weight]*0.5; %*(1-params.t);

function res = my_velocity_components(glob,params);
resx = [5*(1-glob(:,2).^2),...
        zeros(size(glob,1),1)];
%resy = [zeros(1,size(glob,2));...
%        -2*((1-glob(1,:).^2 .*(glob(1,:)>=0) & (glob(1,:)<=1)...
%	     ))];
resy = [zeros(size(glob,1),1),...
        -1*((4-glob(:,1).^2))];
res = {resx, resy};

%function res = my_velocity_coefficients2(params);
%res = 1;

%function res = my_velocity_components2(glob,params);
%res = {repmat([0;-1],1,size(glob,2))};

%| \docupdate 
