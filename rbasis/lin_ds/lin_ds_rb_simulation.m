function rb_sim_data = lin_ds_rb_simulation(model,reduced_data)
%function rb_sim_data = lin_ds_rb_simulation(model,reduced_data)
%
% function performing a reduced basis simulation with a
% theta-scheme time discretization.

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


% Bernard Haasdonk 2.4.2009

if model.verbose >= 9
  disp('entered reduced simulation ');
end;

%keyboard;

% possibly shrink reduced_data detected by this routine 
reduced_data_t = lin_ds_reduced_data_subset(model,reduced_data);

%model.affine_decomp_mode = 'complete';
model.decomp_mode = 2; % coefficients

model = model.set_time(model,0);
%model.t = 0;

if isfield(model,'save_time_indices')
  save_time_index = zeros(1,model.nt+1);
  save_time_index(model.save_time_indices(:)+1) = 1;
  save_time_indices = find(save_time_index==1)-1;
else
  save_time_index = ones(1,model.nt+1);
  save_time_indices = 1:(model.nt+1);
end;

x0coeff = model.x0_function_ptr(model,[]);
xr = lincomb_sequence(reduced_data_t.x0r,x0coeff);


% instead of large dense matrix allocation equivalent:
% a' Mtemp b = 
% a' G b - (a' W) (V' G b) - (a' G V) (W' b) + (a' W) (V' G V) (W' b)
%
% these components can not be summed up here,
% as we must be able to extract online data for other reduced dimensions!!!
% hence no longer computation of m here 
% but decomposition as follows:
%
% m00: 2d cell array of values x0{q}'G x0{q'} 
% VtGx0: 1d cell array of of vectors  (V' G x0{q})
% Wtx0:  1d cell array of of vectors  (W' x0{q})
% VtGV : matrix V' G V

if model.enable_error_estimator
  %e0sqr = lincomb_sequence2(reduced_data_t.m,x0coeff,x0coeff);
  x0tGx0 = lincomb_sequence2(reduced_data_t.m00,x0coeff,x0coeff);
  Wtx0 = lincomb_sequence(reduced_data_t.Wtx0,x0coeff);
  VtGx0 = lincomb_sequence(reduced_data_t.VtGx0,x0coeff);
  VtGV = reduced_data_t.VtGV; 
  e0sqr = x0tGx0 - (Wtx0)'*VtGx0 - (VtGx0)' * Wtx0 + (Wtx0)' * VtGV * Wtx0;
  Rnormsqrs = zeros(1,length(save_time_indices));
end;

D = model.D_function_ptr(model,[]);
Xr = zeros(length(xr(:)),length(save_time_indices));
Y = zeros(size(reduced_data_t.Cr{1},1),size(Xr,2));

if model.enable_error_estimator
  DeltaX = zeros(1,length(save_time_indices));
end;

time_indices = zeros(1,length(save_time_indices));
time =  zeros(1,length(save_time_indices));

t_column = 1; % next column index to be filled in output
t_ind = 0; % current t_ind between 0 and nt
t = 0; % absolute time between 0 and T

if model.enable_error_estimator
  DeltaXt = model.state_bound_constant_C1*sqrt(max(e0sqr,0));
end;

if save_time_index(t_ind+1)
  Xr(:,t_column) = xr(:);
  time_indices(t_column) = t_ind;
  time(t_column) = t;
  %u = model.u_function_ptr(model);  

  Ccoeff = model.C_function_ptr(model,[]);
  Cr = lincomb_sequence(reduced_data_t.Cr,Ccoeff);
  u = model.u_function_ptr(model);  
  Y(:,t_column) = Cr * xr + D * u;

  if model.enable_error_estimator
    DeltaX(t_column) = DeltaXt;
  end;
  t_column = t_column + 1;  
end;

%keyboard;

model.dt = model.T/model.nt;
dt = model.dt;
theta = model.theta;

operators_required = 1;
if ~isfield(model,'data_const_in_time')
  model.data_const_in_time = 0;
end;
previous_implicit_operators_available = 0;

% loop over time
for t_ind = 1:model.nt  
%  keyboard;
  if model.verbose >= 10
    disp(['entered time-loop step ',num2str(t_ind)]);
  end;
  if model.verbose == 9
    fprintf('.');
  end;

  t = (t_ind-1)*dt;
  model = model.set_time(model,t);
%  model.t = t;
  
  if operators_required

    if previous_implicit_operators_available
      Ar0 = Ar1;
      Br0 = Br1;
    else
      A0coeff = model.A_function_ptr(model,[]);
      Ar0 = lincomb_sequence(reduced_data_t.Ar,A0coeff);
      B0coeff = model.B_function_ptr(model,[]);
      Br0 = lincomb_sequence(reduced_data_t.Br,B0coeff);
    end;
    
    if model.enable_error_estimator
      M1 = lincomb_sequence2(reduced_data_t.M1,A0coeff,A0coeff);  
      M2 = lincomb_sequence2(reduced_data_t.M2,B0coeff,B0coeff);  
      M3 = reduced_data_t.M3;
      M4 = lincomb_sequence2(reduced_data_t.M4,B0coeff,A0coeff);  
      M5 = lincomb_sequence(reduced_data_t.M5,A0coeff);  
      M6 = lincomb_sequence(reduced_data_t.M6,B0coeff);       
      % output matrices required at next time!!
    end;
    model = model.set_time(model,t+dt);
    if theta > 0
      A1coeff = model.A_function_ptr(model,[]);
      Ar1 = lincomb_sequence(reduced_data_t.Ar,A1coeff);
      B1coeff = model.B_function_ptr(model,[]);
      Br1 = lincomb_sequence(reduced_data_t.Br,B1coeff);
      previous_implicit_operators_available = 1;
    end;
    %    model.t = t+model.dt;
    C1coeff = model.C_function_ptr(model,[]);
    Cr1 = lincomb_sequence(reduced_data_t.Cr,C1coeff);
    D1 = model.D_function_ptr(model,[]);
    model = model.set_time(model,t);
%    model.t = t;
    if model.data_const_in_time
      operators_required = 0;
    end;
  end;
  
  u0 = model.u_function_ptr(model);  
  xr_old = xr;
  ddt_xr = Ar0 * xr + Br0 * u0;
  if theta == 0
    xr = xr + dt * ddt_xr;
  else
    model = model.set_time(model,t+dt);
    u1 = model.u_function_ptr(model);  
    model = model.set_time(model,t+dt);
    % solve M1 x1= rhs := M0 x0 + b
    Ma1 = eye(size(Ar0))-dt * theta * Ar1;
    Ma0 = eye(size(Ar0))+dt * (1-theta) * Ar0;
    b = dt * ( theta * Br1 * u1  + (1-theta)*Br0 *u0);
    rhs = Ma0 * xr + b;
    xr = Ma1 \ rhs;
  end;
  
  if model.enable_error_estimator
    
    % residual at previous timestep => "rectangular" integration
    Rnormsqr = xr_old' * M1 * xr_old + u0' * M2 * u0 + ...
	ddt_xr' * M3 * ddt_xr + ...
	2*u0' * M4 * xr_old - 2*ddt_xr' * M5 * xr_old - ...
	2 * ddt_xr' * M6 * u0;
    % error estimator at next time
    DeltaXt = DeltaXt + ...
	      model.state_bound_constant_C1 * model.dt ...
	      * sqrt(max(Rnormsqr,0));
%    keyboard;
  end;
    
  % store result
  if save_time_index(t_ind+1)
    Xr(:,t_column) = xr;
    time_indices(t_column) = t_ind;    
    % for output u required at next timestep!
    model = model.set_time(model,t+dt);
%    model.t = t+model.dt;
    u = model.u_function_ptr(model);  
    model = model.set_time(model,t);
%    model.t = t;
    Y(:,t_column) = Cr1 * xr + D1 * u;    
    time(t_column) = t;

    if model.enable_error_estimator
      DeltaX(t_column) = DeltaXt;
      Rnorms(t_column) = sqrt(max(Rnormsqr,0));
    end;
    
    t_column = t_column + 1;  
  end;

end;

if model.verbose == 9
  fprintf('\n');
end;

rb_sim_data = [];
rb_sim_data.Xr = Xr;
rb_sim_data.Y = Y;
if model.enable_error_estimator
  rb_sim_data.DeltaX = DeltaX;
  rb_sim_data.DeltaY = DeltaX * model.output_bound_constant_C2;
  rb_sim_data.Rnorms = Rnorms;
end;
rb_sim_data.time = time;
rb_sim_data.time_indices = time_indices;

%if model.debug
%  if max(rb_sim_data.Y)>1e10
%    error('stability problem, please inspect!');
%  end;
%end;

%keyboard;

%if model.debug
%  disp('short before leavin lin_ds_rb_simulation');
%end;
