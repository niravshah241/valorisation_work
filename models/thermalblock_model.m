function model = thermalblock_model(params)
% Thermal Block example similar as described in the book of 
% A.T. patera and G. Rozza (just one parameter more)
%
% i.e. 
% - div ( a(x) grad u(x)) = q(x)    on Omega
%                   u(x)) = g_D(x)  on Gamma_D
%     a(x) (grad u(x)) n) = g_N(x)  on Gamma_N
%                       s = l(u) linear output functional
%
% where Omega = [0,1]^2 is divided into B1 * B2 blocks 
% QA := B1*B2. The heat conductivities are given by mu_i:
%
%   -------------------------
%   | ...   | ..... | mu_QA |
%   -------------------------
%   | ....  | ...   | ..... |
%   -------------------------
%   | mu_1  | ...   | mu_B1 |
%   -------------------------
%
% `Gamma_D =` upper edge
% `Gamma_N = boundary(Omega)\Gamma_D`
%
% `a(x) = mu_i(x) if x\in B_i`
% `q(x) = 0`
% `g_D  = 0` on top
% `g_N  = 1` on lower edge, 0 otherwise
% `l(u)` = average over lower edge

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


% S. Langhof and B. Haasdonk 26.2.2011
% I. Maier 26.04.2011

model = lin_stat_model_default;
model.name = 'thermalblock';

if nargin == 0
    params.B1=3;
    params.B2=3;
%    model.default = '3x3 demo version';
end
% number of blocks in X and Y direction
model.B1 = params.B1;
model.B2 = params.B2;
model.number_of_blocks = params.B1*params.B2;

% each block has a different heat conductivity
% upper left has conductivity 1
mu_names = {};
mu_ranges = {};
mu_range = [0.1,10];
for p = 1:model.number_of_blocks
  mu_names = [mu_names,{['mu',num2str(p)]}];
  mu_ranges = [mu_ranges,{mu_range}];
end;
model.mu_names = mu_names;
model.mu_ranges = mu_ranges;

%default values 1 everywhere
model.mus = ones(model.number_of_blocks,1);
model.get_mu = @(model) model.mus(:);
model.set_mu = @thermalblock_set_mu;

% simple snapshots without orthogonalization
model.RB_generation_mode = 'lagrangian';
model.RB_mu_list = {[0.1,1,1,1,1,1,1,1,1],...
		    [1,0.1,1,1,1,1,1,1,1],...
		    [1,1,0.1,1,1,1,1,1,1],...
		    [1,1,1,0.1,1,1,1,1,1],...
		    [1,1,1,1,0.1,1,1,1,1],...
		    [1,1,1,1,1,0.1,1,1,1],...
		    [1,1,1,1,1,1,0.1,1,1],...
		    [1,1,1,1,1,1,1,0.1,1],...
		    [1,1,1,1,1,1,1,1,0.1]};

%%%%%%%%% set data functions
model.has_diffusivity = 1;
model.has_output_functional = 1;
model.has_dirichlet_values = 1;
model.has_neumann_values = 1;

% zero dirichlet values, i.e. 1 component, Q_dir = 1
dirichlet_values_coefficients = @(dummy,params) [0];
dirichlet_values_components = @(glob,params) {zeros(size(glob,1),1)};
model.dirichlet_values = @(glob,params) ...
    eval_affine_decomp_general(dirichlet_values_components, ...
			       dirichlet_values_coefficients,glob,params);

% 1/0 neumann values depending on edge, non parametric, i.e Q_neu = 1;
neumann_values_coefficients = @(dummy,params) 1; % single component
model.neumann_values = @(glob,params) ...
    eval_affine_decomp_general(@thermalblock_neumann_values_components, ...
			       neumann_values_coefficients, glob, params);

% diffusion tensor: each row four entries a11,a_21,a_12,a_22. 
% a11(x)=a22(x) = mu_i if x in block i, a12=a21 = 0. 
diffusivity_tensor_coefficients = @(dummy,params) ...
    params.mus(:);
model.diffusivity_tensor = @(glob,params) ...
    eval_affine_decomp_general(...
	@thermalblock_diffusivity_tensor_components, ...
	diffusivity_tensor_coefficients, glob,params);

% only useful for detailed simulation or nonlinear outputs:
%model.output_functional = @output_functional_boundary_integral;

model.operators_output = @fem_operators_output_boundary_integral;
% output weight function: simply 1 on lower edge, 0 else => Ql = 1 component
output_function_components = @(glob,model) {1.0*(glob(:,2)<eps)};
output_function_coefficients = @(glob,model) 1;
model.output_function = @(glob,params) ...
    eval_affine_decomp_general(...
	output_function_components,...
	output_function_coefficients, glob,params);

model.output_integral_qdeg = 2; %

%%%%%%%%%% discretization settings

if ~isfield(params,'numintervals_per_block')
  params.numintervals_per_block = 5;
end;

model.gridtype = 'triagrid';
model.xnumintervals = params.numintervals_per_block*params.B1;
model.ynumintervals = params.numintervals_per_block*params.B2;
model.xrange = [0,1];
model.yrange = [0,1];

if ~isfield(params,'pdeg')
  params.pdeg = 1; 
end;

if ~isfield(params,'qdeg')
  params.qdeg = 2; 
end;

% numerics settings:
model.pdeg = params.pdeg; % degree of polynomial functions
model.qdeg = params.qdeg; % quadrature degree
model.dimrange = 1; % scalar solution

model.boundary_type = @thermalblock_boundary_type;
%model.normals = @my_normals;

% new plot routine, additional block boundaries are plotted
model.plot_sim_data = @thermalblock_plot_sim_data;
model.compute_output_functional = 1;
model.yscale_uicontrols = 0.2;

% for error estimation:
model.coercivity_alpha = @(model) min(model.get_mu(model));
model.continuity_gamma = @(model) max(model.get_mu(model));
model.enable_error_estimator = 1;

% local data functions:
model = elliptic_discrete_model(model);

%%%%%%% auxiliary functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function res = thermalblock_boundary_type(glob,params)
res = zeros(size(glob,1),1);
i  = find(glob(:,1)<=1e-10);
res(i) = -2;
i  = find(glob(:,1)>=1-1e-10);
res(i) = -2;
i  = find(glob(:,2)<=1e-10);
res(i) = -2;
i  = find(glob(:,2)>=1-1e-10);
res(i) = -1;

function model = thermalblock_set_mu(model,mu)
if length(mu)~=model.number_of_blocks
  error('length of mu does not fit to number of blocks!');
end;
model.mus = [mu(:)];

function comp = thermalblock_neumann_values_components(glob,params) 
res = zeros(size(glob,1),1);
i = find(glob(:,2)<eps);
res(i) = 1;
comp = {res};

function res = thermalblock_diffusivity_tensor_components(glob,params)
% for each point in glob find global block number
% xi: range 0 ... B1-1
xi = floor(glob(:,1)*params.B1);
i = find(xi>=params.B1);
if ~isempty(i)
  xi(i) = params.B1-1; 
end;
% xi: range 0 ... B1-1
yi = floor(glob(:,2)*params.B2);
i = find(yi>=params.B2);
if ~isempty(i);
  yi(i) = params.B2-1; 
end;
block_index = yi*params.B1+xi+1;
zeroblock = zeros(size(glob,1),4);
res = cell(1,params.number_of_blocks);
for q = 1:params.number_of_blocks;
  block = zeroblock;
  i = find(block_index==q);
  if ~isempty(i)
    block(i,1) = 1;
    block(i,4) = 1;
  end;
  res{q} = block;
end;
% check!
%block = zeroblock;
%for q=1:params.number_of_blocks;
%  block = block + res{q};
%end;
%disp('check diffusivity matrix!')
%keyboard;

function alpha = thermalblock_coercivity_alpha(model)
alpha = min(model.get_mu(model))


function p = thermalblock_plot_sim_data(model,detailed_data,sim_data, ...
				    plot_params)
if nargin<4
  plot_params = [];
end;
if ~isfield(plot_params,'no_lines')
  plot_params.no_lines = 1;
end;
p = lin_stat_plot_sim_data(model,detailed_data,sim_data,plot_params);
% plot coarse mesh
if ~isfield(plot_params,'plot_blocks')
  plot_params.plot_blocks = 1;
end;
if plot_params.plot_blocks
  X = [0:1/model.B1:1];
  Y = [0:1/model.B2:1];
  l1 = line([X;X],...
	    [zeros(1,model.B1+1);...
	     ones(1,model.B1+1)]);
  set(l1,'color',[0,0,0],'linestyle','-.');
  %keyboard;
  l2 = line([zeros(1,model.B2+1);...
	     ones(1,model.B2+1)],...
	    [Y;Y]);
  set(l2,'color',[0,0,0],'linestyle','-.');
  p = [p(:);l1(:);l2(:)];
end;










