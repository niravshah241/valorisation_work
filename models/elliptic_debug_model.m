function model = elliptic_debug_model(params);
%function model = elliptic_debug_model(params);
% 
% small example of a model, i.e. a structure describing the data
% functions and geometry information of a general elliptic equation consisting 
% of diffusion, convection, reaction equation:
%
% - div ( a(x) grad u(x)) + div (b(x) u(x)) + c(x) u(x) = f(x)    on Omega
%                                                 u(x)) = g_D(x)  on Gamma_D
%                                   a(x) (grad u(x)) n) = g_N(x)  on Gamma_N
%          alpha(x) u(x) +  beta(x) a(x) (grad u(x)) n) = g_R(x)  on Gamma_R
%
%                                   s = l(u) linear output functional
% 
%  Here, we denote the functions as
%                   u: scalar 'solution' (if known, for validation purpose)
%                   f: scalar 'source'
%                   a: tensor valued 'diffusivity_tensor'
%                   b: vector valued 'velocity'
%                   c: scalar 'reaction'
%                 g_D: scalar 'dirichlet_values'
%                 g_N: scalar 'neumann_values'
%                 g_R: scalar 'robin_values'
%               alpha: scalar 'robin_alpha'
%                beta: scalar 'robin_beta'
%
% Each function allows the evaluation in many points
% simultaneuously by
%
%        model.source(glob)
% or     model.source(glob,params)
%
% where glob is a n times 2 matrix of row-wise points. The result
% is a n times 1 vector of resulting values of f.
%
%        model.diffusivity_tensor(glob)
% or     model.diffusivity_tensor(glob,params)
%
% where glob is a n times 2 matrix of row-wise points. The result
% is a n times 4 matrix of resulting values of diffusivity, where
% the values of a are sorted in matlab-order as [a_11, a_21, a_12, a_22];
%
% additionally, the model has a function, which determines, whether
% a point lies on a Dirichlet, Neumann or Robin boundary:
%        
%           model.boundary_type(glob)
%                0 no boundary (inner edge or point)
%               -1 indicates Dirichlet-boundary
%               -2 indicates Neumann-boundary
%               -3 indicates Robin-boundary
%
% Additionally, the normals in a boundary point can be requested by
%
%           model.normals(glob)
% Here, glob are assumed to be boundary points lying ON THE
% INTERIOR of an edge, such that the outer unit normal is well-defined.
%
% The latter 2 methods boundary_type() and normals() need not be
% implemented for grid-based methods, as the normals simply can be
% obtained by the grid. The methods are only required, if using
% meshless methods with data functions that use normals. 
%
% The output functional must be specified 
%           s = model.output_functional(model, model_data, u)
% where u is a dof vector of the discrete solution. model_data is
% assumed to contain the grid 
% 
% possible fields of params:
%     numintervals: the unit square is divided into numintervals
%                   intervals per side. default is 10;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% The data functions given in this model are a parametrized model
% to be used for debugging purpose. The solution for all parameters
% is identical, but the data functions change. 
% The model parameters are mu_a in [-1,1], mu_b in [0,1] and mu_c
% in [0,1]. This model should be used for validating any new
% elliptic solver BY PARAMETER VARIATION.
% by mu_b=mu_c = 0 the advection and reaction terms can be turned
% off. by params.all_dirichlet_boundary = 1, all boundary types can
% be set to dirichlet to turn off Neumann and Robin boundaries.
% 
% The data functions
% are
%
%  a = [1,0;0,1] + mu_a [0, 1/2; 1/2, 0];
%  b = [1/2;1/2] (x-y) mu_b; 
%  c = (1-x)mu_c 
% 
% the exact solution is assumed to be u(x) = sin(k_x x) sin(k_y y) 
%
% The computational domain is the unit interval. Dirichlet
% boundaries assigned right and upper. Neumann assigned left and
% Robin assigned at bottom.
%
% The source term, the boundary values are then specified by the PDE.
%
% An output functional is simply the average over the complete domain

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


% B. Haasdonk 24.1.2011

if nargin == 0 
  params = [];
end;

if ~isfield(params,'numintervals')
  params.numintervals = 10; % 2 x 2 grid! with 8 triangles
end;

if ~isfield(params,'pdeg')
  params.pdeg = 1; 
end;

if ~isfield(params,'qdeg')
  params.qdeg = 2; 
end;

model = [];
model = lin_stat_model_default;

model.kx = 4;
model.ky = 4;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Play with the following parameters for debugging a discretization
%%% initial setting: simple poisson problem:
%model.mu_a = 0;
%model.mu_b = 0;
%model.mu_c = 0;
%model.all_dirichlet_boundary = 1;
%%% more complex volume integrals
%model.mu_a = 1;
%model.mu_b = 1;
%model.mu_c = 1;
%model.all_dirichlet_boundary = 1;
%%% robin and neumann boundary conditions and complex volume integrals
model.mu_a = 1;
model.mu_b = 1;
model.mu_c = 1;
model.all_dirichlet_boundary = 0;

% parameter settings for rb-methods:
model.mu_names ={'mu_a','mu_b','mu_c'};
model.mu_ranges ={[0,1],[0,1],[0,1]};
model.RB_mu_list = {[0,0,0],[0,0,1],[0,1,0],[0,1,1],[1,0,0],...
		    [1,0,1],[1,1,0],[1,1,1]};
model.RB_generation_mode = 'lagrangian';

%%%%%%% data function specification:

% if any of the following flags is set, the corresponding function
% pointers must be provided below.

model.has_reaction = 1;
%model.has_reaction = 0;
model.has_source = 1;
model.has_reaction = 1;
%model.has_reaction = 0;
model.has_advection = 1;
%model.has_advection = 0;
model.has_diffusivity = 1;
model.has_output_functional = 1;
model.has_dirichlet_values = 1;
model.has_neumann_values = 1;
%model.has_neumann_values = 0;
%model.has_robin_values = 1;
model.has_robin_values = 1;

% solution is known: for validation purpose
model.solution = @(glob,params) ...
    sin(glob(:,1)*params.kx).* ...
    sin(glob(:,2)*params.ky);
% gradient: each row one gradient
model.solution_gradient = @(glob,params) ...
    [params.kx*cos(glob(:,1)*params.kx).* sin(glob(:,2)*params.ky), ...
     params.ky*sin(glob(:,1)*params.kx).* cos(glob(:,2)*params.ky) ...
    ];
% hessian: each row: u_xx u_yx, u_xy, u_yy
model.solution_hessian = @(glob,params) ...
    [-params.kx*params.kx*sin(glob(:,1)*params.kx).* sin(glob(:,2)*params.ky), ...
     params.kx * params.ky * cos(glob(:,1)*params.kx).* cos(glob(:,2)*params.ky) ...
     params.kx* params.ky * cos(glob(:,1)*params.kx).* cos(glob(:,2)*params.ky), ...
     -params.ky*params.ky * sin(glob(:,1)*params.kx).* sin(glob(:,2)*params.ky) ...
    ];

% params is assumed to be the model
% c = (1-x) mu_c + x (1-mu_c)
reaction_coefficients = @(dummy,params) ...
    [params.mu_c; 1-params.mu_c];
reaction_components = @(glob,params) ...
    {1-glob(:,1), glob(:,1)};
model.reaction = @(glob,params) ...
    eval_affine_decomp_general(reaction_components, ...
		       reaction_coefficients,glob,params);

% b = [1/2; 1/2] *x * mu_b, single component
velocity_coefficients = @(dummy,params) ...
    [params.mu_b];
velocity_components = @(glob,params) ...
    {ones(size(glob))*0.5.*repmat([glob(:,1)-glob(:,2)],1,2)};
model.velocity = @(glob,params) ...
    eval_affine_decomp_general(velocity_components, ...
		       velocity_coefficients, glob, params);

% a = eye(2) +  mu_a [0,1/2; 1/2, 0], two components
diffusivity_tensor_coefficients = @(dummy,params) ...
    [1; params.mu_a];
diffusivity_tensor_components = @(glob,params) ...
    {...
	[ones(size(glob,1),1),...
	 zeros(size(glob,1),1),...
	 zeros(size(glob,1),1),...
	 ones(size(glob,1),1)],...
	0.5 * [zeros(size(glob,1),1),...
	       ones(size(glob,1),1),...
	       ones(size(glob,1),1),...
	       zeros(size(glob,1),1)] ...
    };
model.diffusivity_tensor = @(glob,params) ...
    eval_affine_decomp_general(diffusivity_tensor_components, ...
		       diffusivity_tensor_coefficients, glob,params);

% dirichlet values:
% simply exact solution, not parameter dependent, 1 component!
dirichlet_values_coefficients = @(dummy,params) [1];
dirichlet_values_components = @(glob,params) ...
    {model.solution(glob,params)};
model.dirichlet_values = @(glob,params) ...
    eval_affine_decomp_general(dirichlet_values_components, ...
		       dirichlet_values_coefficients,glob,params);

% neumann_values: 
% simply exact solution, multiplied with diffusion-tensor and
% normal, hence coefficients and component number coincide
neumann_values_coefficients = @(dummy,params) ...
    diffusivity_tensor_coefficients(dummy,params);
model.neumann_values = @(glob,params) ...
    eval_affine_decomp_general(@my_neumann_values_components, ...
		       neumann_values_coefficients, glob, params);

% robin_values: 
% simply weighted sum of dirichlet and neumann values
model.robin_alpha = @(glob,params) ones(size(glob,1),1);
model.robin_beta = @(glob,params) ones(size(glob,1),1);
% first dirichlet coefficients, then neumann coefficients
model.robin_values = @(glob,params) ... 
    eval_affine_decomp_general(@my_robin_values_components, ...
		       @my_robin_values_coefficients,glob,params);

% source terms: obtain from data functions
model.source = @(glob,params) ...
    eval_affine_decomp_general(@my_source_components, ...
		       @my_source_coefficients,glob,params);

% output functional, e.g. average over unit-square:
model.output_functional = @output_functional_volume_integral;
model.output_function = @output_function_box_mean;
model.sbox_xmin = 0;
model.sbox_ymin = 0;
model.sbox_xmax = 1;
model.sbox_ymax = 1;
model.output_integral_qdeg = 2; %

%%%%%%% geometry specification and discretization:

model.gridtype = 'triagrid';
model.xnumintervals = params.numintervals;
model.ynumintervals = params.numintervals;
model.xrange = [0,1];
model.yrange = [0,1];

% numerics settings:
model.pdeg = params.pdeg; % degree of polynomial functions
model.qdeg = params.qdeg; % quadrature degree
model.dimrange = 1; % scalar solution

% The following 2 methods boundary_type() and normals() need not be
% implemented for grid-based methods, as the normals simply can be
% obtained by the grid. The methods are only required, if using
% meshless methods with data functions that use normals. 
if model.all_dirichlet_boundary == 1
  model.boundary_type = @all_dirichlet_boundary_type;
else
  model.boundary_type = @mixed_boundary_type;
end;
model.normals = @my_normals;

% additional settings for reduced basis approach:
model.mu_names = {'mu_a','mu_b','mu_c'};
model.mu_ranges = {[-1,1],[0,1],[0,1]};

%%%%%%% auxiliary functions:

% all edges of unit square are dirichlet, other inner
function res = all_dirichlet_boundary_type(glob,params)
res = zeros(size(glob,1),1);
i  = find(glob(:,1)<=1e-10);
i  = [i, find(glob(:,1)>=1-1e-10)];
i  = [i, find(glob(:,2)<=1e-10)];
i  = [i, find(glob(:,2)>=1-1e-10)];
res(i) = -1;

% right and upper edges of unit square are dirichlet, 
% left neumann, lower robin
function res = mixed_boundary_type(glob,params)
res = zeros(size(glob,1),1);
i  = find(glob(:,1)<=1e-10);
res(i) = -2;
i  = find(glob(:,2)<=1e-10);
res(i) = -3;
i  = find(glob(:,1)>=1-1e-10);
res(i) = -1;
i  = find(glob(:,2)>=1-1e-10);
res(i) = -1;

function res = my_normals(glob,params);
res = zeros(size(glob,1),2); % each row one normal
i  = find(glob(:,1)>1-1e-10);
res(i,1)= 1.0;
i  = find(glob(:,1)<1e-10);
res(i,1)= -1.0;
i  = find(glob(:,2)>1-1e-10);
res(i,2)= 1.0;
i  = find(glob(:,2)<1e-10);
res(i,2)= -1.0;
% remove diagonal normals
i = find (sum(abs(res),2)>1.5);
res(i,1)= 0;

% neumann boundary value components: 
% multiplication of diffusivity components with solution gradient
function res = my_neumann_values_components(glob,params);
diff_comp = params.diffusivity_tensor(glob,params);
u_grad = params.solution_gradient(glob,params);
normals = params.normals(glob,params);
res = cell(1,length(diff_comp));
for q = 1:length(diff_comp)
  % Agradu: each row one Agradu value
  if ~isempty(diff_comp{q})
    Agradu = [diff_comp{q}(:,1).* u_grad(:,1) + ...
	      diff_comp{q}(:,3).* u_grad(:,2), ... 
	      diff_comp{q}(:,2).* u_grad(:,1) + ...
	      diff_comp{q}(:,4).* u_grad(:,2)]; 
    res{q} = Agradu(:,1).* normals(:,1) + Agradu(:,2).* normals(:,2); 
  else
    res{q} = zeros(size(glob,1),1);
  end;
end;

function res = my_robin_values_coefficients(dummy,params)
res  = [params.dirichlet_values(dummy,params); ...
	params.neumann_values(dummy,params)];

function res = my_robin_values_components(glob,params)
dir_comp = params.dirichlet_values(glob,params);
neu_comp = params.neumann_values(glob,params);
robin_alpha = params.robin_alpha(glob,params);
robin_beta = params.robin_beta(glob,params);
for q = 1:length(dir_comp)
  dir_comp{q} = dir_comp{q}.*robin_alpha;
end;
for q = 1:length(neu_comp)
  neu_comp{q} = neu_comp{q}.*robin_beta;
end;
res =  [dir_comp, neu_comp];

% source term: diffusivity, then advection, then reaction contributions
function res = my_source_coefficients(dummy,params);
res = [-params.diffusivity_tensor(dummy,params);...
       params.velocity(dummy,params);...
       params.reaction(dummy,params)];

% source term: diffusivity, then advection, then reaction
% contribution. Assume that a(x) has vanishing derivative, 
% i.e piecewise constant!!!
% velocity assumed to be divergence free!
function res = my_source_components(glob,params);
u  = params.solution(glob,params);
u_grad = params.solution_gradient(glob,params);
% hessian: each row: u_xx u_yx, u_xy, u_yy
u_hessian = params.solution_hessian(glob,params);
diff_comp = params.diffusivity_tensor(glob,params);
vel_comp = params.velocity(glob,params);
reac_comp = params.reaction(glob,params);
source_diff_comp = cell(1,length(diff_comp));
source_adv_comp = cell(1,length(vel_comp));
source_reac_comp = cell(1,length(reac_comp));
for q = 1:length(reac_comp)
  source_reac_comp{q} = reac_comp{q}.* u; 
end;
for q = 1:length(vel_comp)
  source_adv_comp{q} = vel_comp{q}(:,1).* u_grad(:,1) + ...
      vel_comp{q}(:,2).* u_grad(:,2);
end;
for q = 1:length(diff_comp)
  source_diff_comp{q} = ...
      diff_comp{q}(:,1).* u_hessian(:,1) + ...
      diff_comp{q}(:,2).* u_hessian(:,2) + ...
      diff_comp{q}(:,3).* u_hessian(:,3) + ...
      diff_comp{q}(:,4).* u_hessian(:,4);
end;
res = [source_diff_comp, source_adv_comp, source_reac_comp];
