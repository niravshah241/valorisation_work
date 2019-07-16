function rb_sim_data = lin_ds_rb_simulation_euler_forward(model,reduced_data)
%function rb_sim_data = lin_ds_rb_simulation_euler_forward(model,reduced_data)
%
% simulation of reduced linear system with explicit Euler forward
% discretization. This method is deprecated, as same result can be
% obtained with theta-scheme, theta=0 in lin_ds_rb_simulation

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

if model.error_estimation
  %e0sqr = lincomb_sequence2(reduced_data_t.m,x0coeff,x0coeff);
  x0tGx0 = lincomb_sequence2(reduced_data_t.m00,x0coeff,x0coeff);
  Wtx0 = lincomb_sequence(reduced_data_t.Wtx0,x0coeff);
  VtGx0 = lincomb_sequence(reduced_data_t.VtGx0,x0coeff);
  VtGV = reduced_data_t.VtGV; 
  e0sqr = x0tGx0 - (Wtx0)'*VtGx0 - (VtGx0)' * Wtx0 + (Wtx0)' * VtGV * Wtx0;
end;

D = model.D_function_ptr(model,[]);
Xr = zeros(length(xr(:)),length(save_time_indices));
Y = zeros(size(reduced_data_t.Cr{1},1),size(Xr,2));

if model.error_estimation
  DeltaX = zeros(1,length(save_time_indices));
end;

time_indices = zeros(1,length(save_time_indices));
time =  zeros(1,length(save_time_indices));

t_column = 1; % next column index to be filled in output
t_ind = 0; % current t_ind between 0 and nt
t = 0; % absolute time between 0 and T

if save_time_index(t_ind+1)
  Xr(:,t_column) = xr(:);
  time_indices(t_column) = t_ind;
  time(t_column) = t;
  %u = model.u_function_ptr(model);  

  Ccoeff = model.C_function_ptr(model,[]);
  Cr = lincomb_sequence(reduced_data_t.Cr,Ccoeff);
  u = model.u_function_ptr(model);  
  Y(:,t_column) = Cr * xr + D * u;

  if model.error_estimation
    DeltaXt = model.state_bound_constant_C1*sqrt(max(e0sqr,0));
    DeltaX(t_column) = DeltaXt;
  end;
  t_column = t_column + 1;  
end;

model.dt = model.T/model.nt;

operators_required = 1;
if ~isfield(model,'data_const_in_time')
  model.data_const_in_time = 0;
end,

% loop over time
for t_ind = 1:model.nt  
%  keyboard;
  if model.verbose >= 10
    disp(['entered time-loop step ',num2str(t_ind)]);
  end;
  if model.verbose == 9
    fprintf('.');
  end;

  t = (t_ind-1)*model.dt;
  model = model.set_time(model,t);
%  model.t = t;
  
  if operators_required
    Acoeff = model.A_function_ptr(model,[]);
    Ar = lincomb_sequence(reduced_data_t.Ar,Acoeff);
    Bcoeff = model.B_function_ptr(model,[]);
    Br = lincomb_sequence(reduced_data_t.Br,Bcoeff);

    if model.error_estimation
      M1 = lincomb_sequence2(reduced_data_t.M1,Acoeff,Acoeff);  
      M2 = lincomb_sequence2(reduced_data_t.M2,Bcoeff,Bcoeff);  
      M3 = reduced_data_t.M3;
      M4 = lincomb_sequence2(reduced_data_t.M4,Bcoeff,Acoeff);  
      M5 = lincomb_sequence(reduced_data_t.M5,Acoeff);  
      M6 = lincomb_sequence(reduced_data_t.M6,Bcoeff);       
      % output matrices required at next time!!
    end;
    model = model.set_time(model,t+model.dt);
%    model.t = t+model.dt;
    Ccoeff = model.C_function_ptr(model,[]);
    Cr = lincomb_sequence(reduced_data_t.Cr,Ccoeff);
    D = model.D_function_ptr(model,[]);
    model = model.set_time(model,t);
%    model.t = t;
    if model.data_const_in_time
      operators_required = 0;
    end;
  end;
  
  u = model.u_function_ptr(model);  
  ddt_xr = Ar * xr + Br * u;
  xr_old = xr;
  xr = xr + model.dt * ddt_xr;

  if model.error_estimation
    
    % residual at previous timestep => "rectangular" integration
    Rnormsqr = xr_old' * M1 * xr_old + u' * M2 * u + ...
	ddt_xr' * M3 * ddt_xr + ...
	2*u' * M4 * xr_old - 2*ddt_xr' * M5 * xr_old - ...
	2 * ddt_xr' * M6 * u;
    
    % error estimator at next time
    DeltaXt = DeltaXt + ...
	      model.state_bound_constant_C1 * model.dt ...
	      * sqrt(max(Rnormsqr,0));
  end;
    
  % store result
  if save_time_index(t_ind+1)
    Xr(:,t_column) = xr;
    time_indices(t_column) = t_ind;    
    % for output u required at next timestep!
    model = model.set_time(model,t+model.dt);
%    model.t = t+model.dt;
    u = model.u_function_ptr(model);  
    model = model.set_time(model,t);
%    model.t = t;
    Y(:,t_column) = Cr * xr + D * u;    
    time(t_column) = t;

    if model.error_estimation
      DeltaX(t_column) = DeltaXt;
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
if model.error_estimation
  rb_sim_data.DeltaX = DeltaX;
  rb_sim_data.DeltaY = DeltaX * model.output_bound_constant_C2;
end;
rb_sim_data.time = time;
rb_sim_data.time_indices = time_indices;

%if model.debug
%  disp('short before leavin lin_ds_rb_simulation');
%end;
