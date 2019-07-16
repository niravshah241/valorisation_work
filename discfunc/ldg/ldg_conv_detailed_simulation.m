function sim_data = ldg_conv_detailed_simulation(model,model_data,dummy)
%function sim_data = ldg_conv_detailed_simulation(model,model_data,dummy)
%
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
%
% fields of sim_data: U an array of columnwise degrees of freedom-vectors

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


% Bernard Haasdonk 27.8.2009

disp(' to be adjusted!!  ');
keyboard;

if model.verbose >= 9
  disp(['entered detailed simulation ']);
end;

if isempty(model_data)
  model_data = gen_model_data(model);
end;

model.decomp_mode = 0; % == complete;
model.t = 0;

% initial values by midpoint evaluation
U0 = model.init_values_algorithm(model,model_data);

%Uzero = U0;
%Uzero.dofs = 0;

% preallocate U
U = zeros(length(U0),model.nt+1);
U(:,1)=U0;
model.dt = model.T/model.nt;
%operators_required = 1;
%if ~isfield(model,'data_const_in_time')
%  model.data_const_in_time = 0;
%end,

params = model;
params.axis_equal = 1;
model.plot(U0,model_data.grid,params);
keyboard;

% loop over time-steps
for t = 1:model.nt  
%  keyboard;
  if model.verbose >= 10
    disp(['entered time-loop step ',num2str(t)]);
  end;
  if model.verbose == 9
    fprintf('.');
  end;
  % get matrices and bias-vector 

  if operators_required
    model.t = (t-1)*model.dt;
    [L_I, L_E, b] = model.operators_ptr(model, model_data);
    if model.data_const_in_time
      operators_required = 0;
    end;
  end;
  
  rhs = L_E * U{t}.dofs + b;
  if isequal(L_I, speye(size(L_I)))
    U{t+1}.dofs = rhs;
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
    U{t+1}.dofs = L_I \ rhs;
    
%    if flag>0
%      disp(['error in system solution, solver return flag = ', ...
%	    num2str(flag)]);
%      keyboard;
%    end;
  end;      
end;

if model.verbose == 9
  fprintf('\n');
end;

sim_data.U = U;

%| \docupdate 
