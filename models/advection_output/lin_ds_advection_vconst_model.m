function ds_model  = lin_ds_advection_vconst_model(params)
%function ds_model  = lin_ds_advection_vconst_model(params)
%
% function providing the advection_fv_output_vconst_model with some
% accelerated pointers as a linear dynamical system, i.e. a lin_ds model

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


% Bernard Haasdonk 30.3.2010

lin_evol_model = advection_fv_output_vconst_model(params);
%lin_evol_model.debug = 1;
%lin_evol_model_data = gen_model_data(lin_evol_model);

% set some additional fields in model which will be copied into
% ds model depending on whether constants are known or not:
estimate_constants = 0;
if estimate_constants
  %%% if constant estimaton has to be performed:
  lin_evol_model.estimate_lin_ds_nmu = 4;
  lin_evol_model.estimate_lin_ds_nX = 10;
  lin_evol_model.estimate_lin_ds_nt = 4;
else
  %%% otherwise:
  % the following values are rough bounds, so could be 
  % non-rigorous, then larger search must be performed.
  lin_evol_model.state_bound_constant_C1 = 1;
  lin_evol_model.output_bound_constant_C2 = 1.2916;
end;

ds_model = lin_ds_from_lin_evol_model(lin_evol_model);
% set accelerated data functions from below
ds_model.A_function_ptr = @(model,model_data) ...
    eval_affine_decomp(@A_components,...
		       @A_coefficients,...
		       model,model_data);

ds_model.B_function_ptr = @(model,model_data) ...
    eval_affine_decomp(@B_components,...
		       @B_coefficients,...
		       model,model_data);
%ds_model.u_function_ptr = @(params) 0.1; % define new

%ds_model.theta = 1;
%disp('set theta to 1!!!');
%keyboard;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% auxialiary coefficient functions for acceleration of 
% advection model: explicit implementation of coefficients
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Acomp = A_components(model,model_data) 
Acomp = model_data.A_components;

function Acoeff = A_coefficients(model); 
%model.base_model.decomp_mode = 2;
%model.base_model.t = model.t;
%mu = get_mu(model);
%model.base_model = model.set_mu(model.base_model,mu);
%[L_I_coeff, L_E_coeff, b_coeff] = ...
%    model.base_model.operators_ptr(model.base_model, );
%% L_E = Id + Delta t * A 
%Acoeff = L_E_coeff(2:end)/model.base_model.dt;
%keyboard;
%Acoeff = -model.base_model.dt * ...
%	 model.base_model.velocity_coefficients_ptr(model.base_model);
%if model.t>0
  Acoeff = - [model.vx_weight, model.vy_weight]*0.5; %*(1-model.t);
%  keyboard;
%end;
  
function Bcomp = B_components(model,model_data) 
Bcomp = model_data.B_components;

function Bcoeff = B_coefficients(model); 
% hmmm the coefficients are now called twice in A and B
% should somehow be cached later...
%model.base_model.decomp_mode = 2;
%model.base_model.t = model.t;
%mu = get_mu(model);
%model.base_model = model.set_mu(model.base_model,mu);
%[L_I_coeff, L_E_coeff, b_coeff] = ...
%    model.base_model.operators_ptr(model.base_model, );
%% L_E = Id + Delta t * A 
%Bcoeff2 = b_coeff/model.base_model.dt;
vcoeff = - [model.vx_weight, model.vy_weight]*0.5;%*(1-model.t);
Q_0 = model.base_model.cone_number;
dcoeff = zeros(1,Q_0);
max_pos = model.base_model.cone_weight * (Q_0-1)+1;
t = model.t;
for q = 1:Q_0  
  dcoeff(q) = (1-(max_pos-q))*(1-t) * ((max_pos>=q) && (max_pos < q+1)) ...
      + (1+(max_pos-q))*(1-t) * ((max_pos>=q-1) && (max_pos < q));
end;
v = -(vcoeff(:)*(dcoeff(:)'));
Bcoeff =  v(:)*model.cone_amplitude;
%if model.t>0
%keyboard;
%end;

%function Ccomp = C_components(model,model_data) 
%Ccomp = model_data.C_components;
%
%function Ccoeff = C_coefficients(model); 
%model.base_model.decomp_mode = 2;
%model.base_model.t = model.t;
%mu = get_mu(model);
%model.base_model = model.set_mu(model.base_model,mu);
%Ccoeff = ...
%    model.base_model.operators_output(model.base_model);

