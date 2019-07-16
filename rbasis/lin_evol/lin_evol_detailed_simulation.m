function sim_data = lin_evol_detailed_simulation(model,model_data)
%function sim_data = lin_evol_detailed_simulation(model,model_data)

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


% required fields of model:
% T             : final time 
% nt            : number of time-intervals until T, i.e. nt+1
%                 solution slices are computed
% init_values_algorithm: name of function for computing the
%                 initvalues-DOF with arguments (grid, params)
%                 example: init_values_cog
% operators_algorithm: name of function for computing the
%                 L_E,L_I,b-operators with arguments (grid, params)
%                 example operators_conv_diff
% data_const_in_time : if this optional field is 1, the time
%                 evolution is performed with constant operators, 
%                 i.e. only the initial-time-operators are computed
%                 and used throughout the time simulation.
% compute_output_functional: flag indicating, whether output
%                 functional is to be computed
%
% return fields of sim_data:
%    U : sequence of DOF vectors
%    y : sequence of output functional values (if activated)
%
% optional fields of model:
%   starting_time_step: starting time step for simulation (for example in
%                   t-partition). Default value is 0.
%   stopping_tim_step: stopping time step for simulation
%
% Bernard Haasdonk 27.8.2009

if model.verbose >= 9
  disp('entered detailed simulation ');
end;

if isempty(model_data)
  model_data = gen_model_data(model);
end;

model.decomp_mode = 0; % == complete;

model.dt = model.T/model.nt;
%Fixing values for t-partition
if ~isfield(model,'starting_time_step')
    model.t = 0;
    t_ind_start =0;
    t_ind_stop = model.nt;
else
    t_ind_start = model.starting_time_step;
    t_ind_stop = model.stopping_time_step;
    model.t = (t_ind_start)*model.dt;
    %model.nt = t_ind_stop-t_ind_start+1; %+1?
end



if isfield(model,'save_time_indices')
  valid_save_time_indices1 = model.save_time_indices>=t_ind_start;
  valid_save_time_indices2 = model.save_time_indices<=t_ind_stop;
  valid_save_time_indices = valid_save_time_indices1.*valid_save_time_indices2;
  valid_save_time_indices = find(valid_save_time_indices);
  save_time_index = zeros(1,model.nt+1);
  save_time_index(model.save_time_indices(valid_save_time_indices)+1) = 1;
  save_time_indices = find(save_time_index==1)-1;
else
  save_time_index = ones(1,model.nt+1);
  save_time_indices = 1:(model.nt+1);
end;

model.nt = t_ind_stop-t_ind_start+1;

% initial values by midpoint evaluation
Ut = model.init_values_algorithm(model,model_data);

% preallocate U
U = zeros(length(Ut),length(save_time_indices));
time_indices = zeros(1,length(save_time_indices));
time =  zeros(1,length(save_time_indices));

operators_required = 1;
%if ~isfield(model,'data_const_in_time')
%  model.data_const_in_time = 0;
%end,

% get operator components
if model.affinely_decomposed
  old_decomp_mode = model.decomp_mode;
  model.decomp_mode = 1;
  [L_I_comp, L_E_comp, b_comp] = model.operators_ptr(model, model_data);
  model.decomp_mode = old_decomp_mode;
end;

t_column = 1; % next column index to be filled in output
t_ind = t_ind_start; % t_ind between 0 and nt
t = model.t; % absolute time between 0 and T

% store init data
if save_time_index(t_ind+1)
  U(:,t_column) = Ut;
  time_indices(t_column) = t_ind;
  time(t_column) = t;
  t_column = t_column + 1;  
end;

if ~isfield(model,'compute_output_functional')
  model.compute_output_functional = 0;
end,

%model.plot(U0,model_data.grid,params);

% loop over fv-steps
for t_ind = (t_ind_start):(t_ind_stop-1)
  t = (t_ind)*model.dt;
  if model.verbose >= 10
    disp(['entered time-loop step ',num2str(t_ind)]);
  end;
  if model.verbose == 9
    fprintf('.');
  end;
  % get matrices and bias-vector 
  
  if operators_required
    
    model.t = t;
    
    % new assembly in every iteration
    if ~model.affinely_decomposed
      [L_I, L_E, b] = model.operators_ptr(model, model_data);
    else
      % or simple assembly based on affine parameter decomposition
      % => turns out to be slower for few components
      
      % test affine decomposition
      if model.debug
        disp('test affine decomposition of operators')
        test_affine_decomp(model.operators_ptr,3,1,model,model_data);
        disp('OK')
      end;
      model.decomp_mode = 2;
      [L_I_coeff, L_E_coeff, b_coeff] = model.operators_ptr(model, ...
      						  model_data);
      model.decomp_mode = old_decomp_mode;
      L_I = lincomb_sequence(L_I_comp, L_I_coeff);
      L_E = lincomb_sequence(L_E_comp, L_E_coeff);    
      b = lincomb_sequence(b_comp, b_coeff);
    end;
    
    if model.data_const_in_time
      operators_required = 0;
    end;

  end;

  rhs = L_E * Ut + b;
  if isequal(L_I, speye(size(L_I)))
    Ut = rhs;
  else
    % solve linear system
    %    disp('check symmetry and choose solver accordingly!');
    %    keyboard;
    % nonsymmetric solvers:
    %  [U(:,t+1), flag] = bicgstab(L_I,rhs,[],1000);
    %  [U(:,t+1), flag] = cgs(L_I,rhs,[],1000);
    %
    % symmetric solver, non pd:
    % omit symmlq:
    % see bug_symmlq.mat for a very strange bug: cannot solve identity system!
    % reported to matlab central, but no solution up to now.
    % [U(:,t+1), flag] = symmlq(L_I,rhs,[],1000);
    %
    %  [U(:,t+1), flag] = minres(L_I,rhs,[],1000);
    % symmetric solver, pd:
    %[U(:,t+1), flag] = pcg(L_I,rhs,[],1000);
    % bicgstab works also quite well:
    %[U(:,t+1), flag] = bicgstab(L_I,rhs,[],1000);
    Ut = L_I \ rhs;

    %    if flag>0
    %      disp(['error in system solution, solver return flag = ', ...
    %	    num2str(flag)]);
    %      keyboard;
    %    end;

  end;
  
  if save_time_index(t_ind+2)
    U(:,t_column) = Ut;
    time_indices(t_column) = t_ind;
    time(t_column) = t;
    t_column = t_column + 1;  
  end;
  
end;

if model.verbose == 9
  fprintf('\n');
end;

sim_data.U = U;
sim_data.time = time;
sim_data.time_indices = time_indices;

if model.compute_output_functional
  % get operators
  v = model.operators_output(model,model_data);
  sim_data.y = (v(:)') * U;
end;

%| \docupdate 
