function lin_ds_model = lin_ds_from_lin_evol_model(lin_evol_model) 
%function lin_ds_model = lin_ds_from_lin_evol_model(lin_evol_model) 
%
% construction of lin_ds_model from lin_evol_model,
% assuming that implicit matrix is not parameter dependent and 
% can be multiplied inversely to right-hand side matrix (==Id)

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


% Bernard Haasdonk 7.9.2009

% assuming lin_evol_model as base_model:

model = lin_ds_model_default;

model.base_model = lin_evol_model; % so access to these settings possible

copyfields = ...
    {'data_const_in_time','verbose','T','nt',...
    'mu_names','mu_ranges','debug'};

for i = 1:length(copyfields)
  model = setfield(model,copyfields{i},getfield(lin_evol_model,copyfields{i}));
end;

if_exists_copyfields = ...
    {'save_time_indices',...
     'state_bound_constant_C1',...
     'output_bound_constant_C2',...
     'estimate_lin_ds_nmu',...
     'estimate_lin_ds_nX',....
     'estimate_lin_ds_nt',...
     'error_estimation',...
     'RB_generation_mode',...
     'RB_num_intervals'...
    };

for i = 1:length(if_exists_copyfields)
  if isfield(lin_evol_model,if_exists_copyfields{i})
    model = setfield(model,if_exists_copyfields{i},...
			   getfield(lin_evol_model, ...
				    if_exists_copyfields{i}));
  end;
end;

model.A_function_ptr = @(model,model_data) ...
    eval_affine_decomp(@A_components,...
		       @A_coefficients,...
		       model,model_data);

model.B_function_ptr = @(model,model_data) ...
    eval_affine_decomp(@B_components,...
		       @B_coefficients,...
		       model,model_data);

model.C_function_ptr = @(model,model_data) ...
    eval_affine_decomp(@C_components,...
		       @C_coefficients,...
		       model,model_data);

model.x0_function_ptr = @(model,model_data) ...
    eval_affine_decomp(@x0_components,...
		       @x0_coefficients,...
		       model,model_data);

model.D_function_ptr = @D_func;
model.u_function_ptr = @(model) 1;

model.set_mu = @set_mu_with_base_model;
mu = lin_evol_model.get_mu(lin_evol_model);
model = model.set_mu(model,mu);
model.set_time = @set_time_with_base_model;

% special gen_model_data routine
model.gen_model_data = @lin_ds_from_lin_evol_gen_model_data;

% special plot routine
model.plot_sim_data = @lin_ds_from_lin_evol_plot_sim_data;

% determin model data, i.e. matrix components, matrix G
model.affinely_decomposed = 1;

model.dim_u = 1;
model.dim_y = 1;

%model_data = gen_model_data(model); % store matrix components
%model.dim_x = size(model_data.A_components{1},1);
model.dim_x = lin_evol_model.dim_U;

%model.get_estimator_from_sim_data = @(sim_data) sim_data.DeltaX(end);
% ds_model.get_estimator_from_sim_data = @(sim_data) sim_data.DeltaY(end);

% be sure to determine correct constants once, if anything in the
% problem setting has changed!!

if ~(isfield(model,'state_bound_constant_C1') & ...
     isfield(model,'output_bound_constant_C2'))
  disp('estimating coefficients...');
  [C1,C2] = lin_ds_estimate_bound_constants(model,model_data);
  disp('...finished');
  model.state_bound_constant_C1 = C1;
  model.output_bound_constant_C2 = C2;
end;

lin_ds_model = model;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Acomp = A_components(model,model_data) 
Acomp = model_data.A_components;

function Acoeff = A_coefficients(model); 
model.base_model.decomp_mode = 2;
model.base_model.t = model.t;
%mu = get_mu(model);
%model.base_model = model.base_model.set_mu(model.base_model,mu);
[L_I_coeff, L_E_coeff, b_coeff] = ...
    model.base_model.operators_ptr(model.base_model, []);
% L_E = Id + Delta t * A 
Acoeff = L_E_coeff(2:end)/model.base_model.dt;
%keyboard;

function Bcomp = B_components(model,model_data) 
Bcomp = model_data.B_components;

function Bcoeff = B_coefficients(model); 
% hmmm the coefficients are now called twice in A and B
% should somehow be cached later...
model.base_model.decomp_mode = 2;
% base_model.t should be correctly set by set_time of model...
%model.base_model.t = model.t;
%mu = get_mu(model);
%model.base_model = model.base_model.set_mu(model.base_model,mu);
[L_I_coeff, L_E_coeff, b_coeff] = ...
    model.base_model.operators_ptr(model.base_model, []);
% L_E = Id + Delta t * A 
Bcoeff = b_coeff/model.base_model.dt;

function Ccomp = C_components(model,model_data) 
Ccomp = model_data.C_components;

function Ccoeff = C_coefficients(model); 
model.base_model.decomp_mode = 2;
% base_model.t should be correctly set by set_time of model...
%model.base_model.t = model.t;
%mu = get_mu(model);
%model.base_model = model.base_model.set_mu(model.base_model,mu);
Ccoeff = ...
    model.base_model.operators_output(model.base_model,[]);

function res = D_func(model,model_data); 
res = zeros(model.dim_y, model.dim_u); 

function x0comp = x0_components(model,model_data) 
x0comp = model_data.x0_components;

function res = x0_coefficients(model); 
model.base_model.decomp_mode = 2;
% base_model.t should be correctly set by set_time of model...
%model.base_model.t = model.t;
%mu = get_mu(model);
%model.base_model = model.base_model.set_mu(model.base_model,mu);
res = model.base_model.init_values_algorithm(model.base_model,[]);

function model = set_mu_with_base_model(model,mu)
model.base_model = model.base_model.set_mu(model.base_model,mu);

model = set_mu_default(model,mu);

function model = set_time_with_base_model(model,time)
model.base_model = model.base_model.set_time(model.base_model,time);
model = set_time_default(model,time);









