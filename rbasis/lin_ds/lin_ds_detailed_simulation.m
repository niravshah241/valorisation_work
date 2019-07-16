function sim_data = lin_ds_detailed_simulation(model,model_data)
%function sim_data = lin_ds_detailed_simulation(model,model_data)
%
% function performing a general time evolution of a linear dynamical
% system 
%
%  @f[       \frac{d}{dt} x = A(t,\mu) x + B(t,\mu) u                     @f]
%  @f[          y = C(t,\mu) x + D(t,\mu) u  \quad  \mbox{ for }
%  \quad t \in [0,T]    @f]
%          
%  @f[             x(0) = x_0(\mu)                @f]
%
% by theta-scheme time discretization.
%
%
% required fields of model:
% T             : final time 
% nt            : number of time-intervals until T, i.e. nt+1
%                 solution slices are computed
% theta         : real value in [0,1]
%             -   0   = explicit
%             -   0.5 = Crank-Nicolson
%             -   1   = Implicit discretization
% x0_function_ptr: pointer to function for computing the
%                 initvalues-DOF with arguments (model)
% A_function_ptr: pointer to function for computing the matrix A
% B_function_ptr: pointer to function for computing the matrix B
% C_function_ptr: pointer to function for computing the matrix C
% D_function_ptr: pointer to function for computing the matrix D
% u_function_ptr: pointer to input function u(t)
% data_const_in_time : if this optional field is 1, the time
%                 evolution is performed with constant operators, 
%                 i.e. only the initial-time-operators are computed
%                 and used throughout the time simulation.
% affinely_decomposed: if this optional field is set, model_data is
%           assumed to contain the components of the matrices and online
%           assembly is performed by coefficient-component linear combinaton
% save_time_indices: list of integer time indices between 0 (!) and nt,
%           which are to be computed, stored and returned
% verbose: integer verbosity level
%
% Required fields of model_data:
% A_components: cell array of affine decomposed A-components  
% B_components: cell array of affine decomposed B-components  
% C_components: cell array of affine decomposed C-components  
% x0_components: cell array of affine decomposed x0-components  
% Required fields of model_data:
%
% Generated fields of sim_data:
% X: matrix containing state vectors of solution trajectory as
%    columns
% Y: matrix containing output vectors of trajectory as columns
% time: vector of time-values in [0,T] for columns in X,Y
% time_indices: integer time step numbers for columns in X,Y 
%
% Return values:
% sim_data: structure with simulation data results

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
  disp('entered detailed simulation ');
end;

%model.affine_decomp_mode = 'complete';
model.decomp_mode = 0;
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

% initial values by midpoint evaluation

if ~model.affinely_decomposed
  xt = model.x0_function_ptr(model,model_data);
  C0 = model.C_function_ptr(model,model_data);
else
  old_decomp_mode = model.decomp_mode;
  model.decomp_mode = 2;
  x0_coeff = model.x0_function_ptr(model,[]);
  xt = lincomb_sequence(model_data.x0_components,x0_coeff);
  C0_coeff = model.C_function_ptr(model,[]);
  C0 = lincomb_sequence(model_data.C_components,C0_coeff);
  model.decomp_mode = old_decomp_mode;
end;
D = model.D_function_ptr(model,model_data); % no decomposition required
X = zeros(length(xt(:)),length(save_time_indices));
Y = zeros(size(C0,1),length(save_time_indices));
time_indices = zeros(1,length(save_time_indices));
time =  zeros(1,length(save_time_indices));

t_column = 1; % next column index to be filled in output
t_ind = 0; % current t_ind between 0 and nt
t = 0; % absolute time between 0 and T

% store init data
if save_time_index(t_ind+1)
%  disp('\n OK, saving initial time')
%    keyboard;
  X(:,t_column) = xt;
  time_indices(t_column) = t_ind;
  time(t_column) = t;
  u = model.u_function_ptr(model);  
  Y(:,t_column) = C0 * xt + D * u;
  t_column = t_column + 1;  
end;

model.dt = model.T/model.nt;
dt = model.dt;
theta = model.theta;
operators_required = 1;
previous_implicit_operators_available = 0;

% loop over time-steps
for t_ind = 1:model.nt  
  t = model.dt * (t_ind-1); % compute new x and y at time t
  if model.verbose >= 10
    disp(['entered time-loop step ',num2str(t_ind)]);
  end;
  if model.verbose >= 9
    fprintf('.');
  end;
  % get matrices and bias-vector 
  model = model.set_time(model,t);
%  model.t = t;
  if operators_required
    if ~model.affinely_decomposed
      if previous_implicit_operators_available
	A0 = A1;
	B0 = B1;
      else
	A0 = model.A_function_ptr(model,model_data);
	B0 = model.B_function_ptr(model,model_data);
      end;
	% C and D required at time t+dt;
      model = model.set_time(model,t+model.dt);
%      model.t = t+model.dt;
      C1 = model.C_function_ptr(model,model_data);
      D1 = model.D_function_ptr(model,model_data);
      if theta>0
	A1 = model.A_function_ptr(model,model_data);
	B1 = model.B_function_ptr(model,model_data);
	previous_implicit_operators_available = 1;
      end;
      model = model.set_time(model,t);
      %      model.t = t;
    else
      model.decomp_mode = 2;
      if ~previous_implicit_operators_available
	A0_coeff = model.A_function_ptr(model,[]);
	B0_coeff = model.B_function_ptr(model,[]);
      end;
      model = model.set_time(model,t+model.dt);
      if theta>0
	A1_coeff = model.A_function_ptr(model,[]);
	B1_coeff = model.B_function_ptr(model,[]);
      end;
      %model.t = t+model.dt;
      C1_coeff = model.C_function_ptr(model,[]);
      model.decomp_mode = old_decomp_mode;
      D1 = model.D_function_ptr(model,[]); % no decomposition of D
      model = model.set_time(model,t);
      %model.t = t;
      if previous_implicit_operators_available
	A0 = A1;
	B0 = B1;
      else	  
	A0 = lincomb_sequence(model_data.A_components,A0_coeff);
	B0 = lincomb_sequence(model_data.B_components,B0_coeff);
      end;
      if theta>0
	A1 = lincomb_sequence(model_data.A_components,A1_coeff);
	B1 = lincomb_sequence(model_data.B_components,B1_coeff);
      end;
      C1 = lincomb_sequence(model_data.C_components,C1_coeff);
    end;
    
    if model.data_const_in_time
      operators_required = 0;
    end;
  end;
  u0 = model.u_function_ptr(model);
  if theta == 0
    xt = xt + dt * A0 * xt + dt * B0 * u0;
  else
    % get input at next time
    model = model.set_time(model,t+model.dt);
    u1 = model.u_function_ptr(model);
    model = model.set_time(model,t);
    % solve M1 x1= rhs := M0 x0 + b
    M1 = speye(size(A0))-dt * theta * A1;
    M0 = speye(size(A0))+dt * (1-theta) * A0;
    b = dt * ( theta * B1 * u1  + (1-theta)*B0 *u0);
    rhs = M0 * xt + b;
%    keyboard;
    xt = M1 \ rhs;
  end;
  %  X(:,t_ind+1) = X(:,t_ind) + model.dt * (A * X(:,t_ind) + B * u);
  %Y(:,t_ind+1) = C * X(:,t_ind) + D * u;

  if save_time_index(t_ind+1)
%    disp('\n OK, saving next time');
%    keyboard;
    X(:,t_column) = xt;
    time_indices(t_column) = t_ind;
    % y evaluation requires + dt !!
    model = model.set_time(model,t+model.dt);
    %model.t = t + model.dt;
    u = model.u_function_ptr(model);
    Y(:,t_column) = C1 * xt + D1 * u;
    time(t_column) = t;
    t_column = t_column + 1;  
    fprintf('.');
  end;
  
end;

if model.verbose >= 9
  fprintf('\n');
end;

sim_data = [];
sim_data.X = X;
sim_data.Y = Y;
sim_data.time = time;
sim_data.time_indices = time_indices;
