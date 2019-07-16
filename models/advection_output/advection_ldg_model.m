function model = advection_ldg_model(params);
%function model = advection_ldg_model(params);
% 
% model for the new advection model with functional output
% two overlapping velocity fields and boundary value modification
% are the parameters. Pure explicit discretization of convection with ldg
% discretization 
%
% example of a model not as class but as structure plus function pointers
% goal: library works with structure- or class-models
%
% (Demonstration of base LDG functionality of rbmatlab library)

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

% specification of the time information
model.T             = 1;
model.nt            = 100;
disp('nt to be adjusted!');
model.dt = model.T/model.nt;
model.verbose = 9;

% grid information
%model.grid_initfile = 'rectangle_triagrid.mat';
model.gridtype = 'triagrid';
model.xnumintervals = 40;
model.ynumintervals = 20;
model.xrange = [0,2];
model.yrange = [0,1];
% set completely dirichlet boundary
model.bnd_rect_corner1 = [0,0]-eps;
model.bnd_rect_corner2 = [2,1]+eps;
model.bnd_rect_index   = [-1];

% Main global function pointers expected in every model
model.detailed_simulation = @ldg_conv_detailed_simulation;
model.gen_model_data =      @lin_evol_gen_model_data;

% model can only perform detailed simulations
%model.gen_detailed_data =   @lin_evol_gen_detailed_data;
%model.gen_reduced_data =     @lin_evol_gen_reduced_data;
%model.rb_simulation =       @lin_evol_rb_simulation;
%model.rb_reconstruction =   @lin_evol_rb_reconstruction;

% Dirichlet and Initial data
model.dirichlet_values_ptr = @dirichlet_values_affine_decomposed;
model.dirichlet_values_coefficients_ptr = @my_dirichlet_values_coefficients;
model.dirichlet_values_components_ptr = @my_dirichlet_values_components;
model.init_values_ptr = @init_values_affine_decomposed;
model.init_values_coefficients_ptr = @my_dirichlet_values_coefficients;
model.init_values_components_ptr = @my_dirichlet_values_components;
model.cone_number = 5; % number of components = cones 
model.cone_weight = 1; % leftmost cone with full weight

model.dimrange = 1; % scalar equation
model.pdeg = 2; % polynomial degree of solution and init ansatz function
model.init_values_qdeg = 4; % quadrature degree
model.init_values_algorithm = @disc_init_values;
model.l2project = @ldg_l2project;
model.evaluate = @ldg_evaluate;
model.evaluate_basis = @ldg_evaluate_basis;
model.plot = @ldg_plot;
model.ndofs_per_element = ldg_ndofs_per_element(model);

model.velocity_ptr = @velocity_affine_decomposed;
model.velocity_coefficients_ptr = @my_velocity_coefficients;
model.velocity_components_ptr = @my_velocity_components;
model.convflux = 'advection';

return; % temporary end





model.operators_ptr = @ldg_operators_explicit;
model.L_I_inv_norm_bound = 1;
model.L_E_norm_bound     = 1;

model.flux_linear = 1;
model.vx_weight = 0.5;
model.vy_weight = 0.5;

% RB information:
model.mu_names           = {'cone_weight','vx_weight','vy_weight'};
model.mu_ranges          = {[0 1],[0 1],[0,1]};

%model.neumann_values_ptr = 0;

%model.name_init_values = 'decomp_function_ptr';
% implement below!!
%model.init_values_coefficients_ptr = @my_init_values_coefficients;
%model.init_values_components_ptr = @my_init_values_components;


model.name_diffusive_num_flux = 'none';

%name_diffusive_num_flux = 'gradient';
    % precomputed, divcleaned velocity field
    name_flux               = 'gdl2';                

    lambda = 2.0e-7;   % v = - lambda * grad p
    name_convective_num_flux = 'lax-friedrichs';
    inner_product_matrix_algorithm = @fv_inner_product_matrix;
    verbose = 5;
    
model.rb_init_values        = @rb_init_values;

    
% further method pointers, which are specific to model:
model.orthonormalize             = @model_orthonormalize_gram_schmidt;
model.inner_product              = @fv_inner_product;
model.PCA                        = @model_PCA_fixspace;
model.cached_detailed_simulation = @cached_detailed_simulation;
model.save_detailed_simulations  = @save_detailed_simulations;
model.load_detailed_simulation   = @load_detailed_simulation;
model.use_velocity_matrixfile = 1;
model.divclean_mode = 'none'; % file is already cleaned by optimization

% set matrix-file name for optional generation or reading of file 
model.velocity_matrixfile = ['vel_',    model.name_flux,'_',...
		    num2str(    model.xnumintervals),'x',...
		    num2str(    model.ynumintervals),...
		    '_l',num2str(    model.lambda),'.mat'];

model.lxf_lambda = 1.0194e+003;
model.data_const_in_time = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% auxiliary functions used as pointers above:
%%%%% check later, if they are of general interest to be exported
%%%%% as standalone functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% new function syntax: global coordinates as rows in glob_coord,
% time and parameter variables in params

% dirichlet-data: convex combination of equidistant cones, decaying
% in time
function res = my_dirichlet_values_coefficients(params);
% res is a vector of coefficients
Q_0 = params.cone_number;
res = zeros(1,Q_0)
max_pos = params.cone_weight * (Q_0-1)+1;
t = params.t;
for q = 1:Q_0  
  res(q) = (1-(max_pos-q))*(1-t) * ((max_pos>=q) && (max_pos < q+1)) ...
	   + (1+(max_pos-q))*(1-t) * ((max_pos>=q-1) && (max_pos < q));
end;

function res = my_dirichlet_values_components(glob,params);
% res is a cell-array of row-vectors of global evaluations
Q_0 = params.cone_number;
res = cell(1,Q_0);
delta_cone = 1/(Q_0+1);
cone_pos_x = delta_cone * (1:Q_0);
for q = 1:Q_0  
  res{q} = 1-min(sqrt((glob(:,1)-cone_pos_x(q)).^2+...
	    (glob(:,2)-1).^2)*(Q_0+1),1);
end;

% velocity field: overlap of y- and x parabolic profiles
function res = my_velocity_coefficients(params);
res = [params.vx_weight, params.vy_weight]*(1-params.t);

function res = my_velocity_components(glob,params);
resx = [(1-glob(:,2).^2);...
        zeros(1,size(glob,1))];
resy = [zeros(1,size(glob,1));...
        (1-glob(:,1).^2.*(glob(:,1)>=0 && glob(:,1)<=1))];
res = {resx, resy};
%| \docupdate 
