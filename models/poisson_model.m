function model = poisson_model(params);
%function model = poisson_model(params)
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
% The data functions given in this model are the simple poisson
% equation with Gamma_D = boundary(Omega), Gamma_N = {}, Gamma_R = {}
%
%    -div (grad u) = f
%            u = 0   on   Gamma_D
%
% with exact solution u(x) = 16 x_1(1-x_1)x_2(1-x_2), 
% i.e. f(x) = 32 (x_1 + x_2 - x_1^2 - x_2^2)

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


% B. Haasdonk 23.11.2010

if nargin == 0 
  params = [];
end;

if ~isfield(params,'numintervals')
  params.numintervals = 10; % 2 x 2 grid! with 8 triangles
end;

if ~isfield(params,'pdeg')
  params.pdeg = 2; % 2 x 2 grid! with 8 triangles
end;

if ~isfield(params,'qdeg')
  params.qdeg = 2; % 2 x 2 grid! with 8 triangles
end;

model = [];
model.rb_problem_type = 'lin_elliptic';

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

% solution is known:
model.solution = @(glob,params) ...
            16 * glob(:,1).*(1-glob(:,1)).*glob(:,2).*(1-glob(:,2));
% for debugging: f = 1; u = -1/4 (x^2+y^2)
%model.solution = @(glob,params) ...
%            -0.25 * (glob(:,1).^2+glob(:,2).^2);

% params is an optional parameter, perhaps useful later
model.source = @(glob,params) ...
            32 * (glob(:,1)+ glob(:,2)- glob(:,1).^2 - glob(:,2).^2);
% for debugging: f = 1; u = -0.25* (x^2+y^2)
%model.source = @(glob,params) ...
%    ones(size(glob,1),1);

model.reaction = @(glob,params) zeros(size(glob,1),1);
model.velocity = @(glob,params) zeros(size(glob,1),2);
model.diffusivity_tensor = @(glob,params) ...
    [ones(size(glob,1),1),...
     zeros(size(glob,1),1),...
     zeros(size(glob,1),1),...
     ones(size(glob,1),1)];
%model.diffusivity_gradient = @(glob,params) zeros(size(glob,1),2);
model.reaction = @(glob,params) zeros(size(glob,1),1);

%  Dirichlet everywhere, see function below.
%model.dirichlet_values = @(glob,params) zeros(size(glob,1),1);
model.dirichlet_values = @(glob,params) ...
    params.solution(glob,params);
model.neumann_values = @(glob,params) zeros(size(glob,1),1);
model.robin_values = @(glob,params) zeros(size(glob,1),1);
model.robin_alpha = @(glob,params) zeros(size(glob,1),1);
model.robin_beta = @(glob,params) zeros(size(glob,1),1);

% output functional, e.g. average over unit-square:
model.output_functional = @output_functional_volume_integral;
model.output_function = @output_function_box_mean;
model.sbox_xmin = 0;
model.sbox_ymin = 0;
model.sbox_xmax = 1;
model.sbox_ymax = 1;
model.output_integral_qdeg = 0; % midpoint integration
% later, e.g. for thermal block:
%model.output_functional = @output_functional_boundary_integral;
%model.output_function = @output_function_dirichlet_indicator;

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
model.boundary_type = @my_boundary_type;
model.normals = @my_normals;

%%%%%%% auxiliary functions:

% all edges of unit square are dirichlet, other inner
% note this function is not used by grid-based methods
function res = my_boundary_type(glob,params)
res = zeros(size(glob,1),1);
i  = find(glob(:,1)<=1e-10);
i  = [i, find(glob(:,1)>=1-1e-10)];
i  = [i, find(glob(:,2)<=1e-10)];
i  = [i, find(glob(:,2)>=1-1e-10)];
res(i) = -1;

function res = my_normals(glob,params);
res = zeros(size(glob,1),2); % each row one normal
i  = find(glob(:,1)>1-1e-10);
res(i,1)= 1.0;
i  = find(glob(:,1)<-1+1e-10);
res(i,1)= -1.0;
i  = find(glob(:,2)>1-1e-10);
res(i,2)= 1.0;
i  = find(glob(:,2)<-1+1e-10);
res(i,2)= -1.0;
% remove diagonal normals
i = find (sum(abs(res),2)>1.5);
res(i,1)= 0;



