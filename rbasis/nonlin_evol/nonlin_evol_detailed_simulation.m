function [sim_data, fl] = nonlin_evol_detailed_simulation(model,model_data)
%function sim_data = nonlin_evol_detailed_simulation(model,model_data)
% computes the detailed numerical scheme as described by 'model' for the
% parameter set in the same data structure 'model'
%
% Return values:
%  sim_data: structure holding the generated Dof Vector
%  fl:       debugging output (second return value from 'L_E_local_ptr')
%
% Generated fields of sim_data:
%  U:        Dof Vector of result for each time instance (matrix of size
%            'H x (nt+1)')
%
% required fields of model:
% T             : final time
% nt            : number of time-intervals until 'T', i.e. 'nt+1' solution
%                 slices are computed
% init_values_algorithm : name of function for computing the initvalues-DOF
%                 with arguments '(model, model_data)' example: fv_init_values()
% L_E_local_ptr :  function pointer to the local space-discretization operator
%                  evaluation with syntax
%                  @code
%                   INC_local = L_local(U_local_ext, ind_local,
%                                      grid_local_ext)
%                  @endcode
%                  where '*_local indicates' results on a set of elements,
%                  '*_local_ext' includes information including their
%                  neighbours, i.e. the grid must be known also for the
%                  neighbours, the values of the previous timstep must be known
%                  on the neighbours. The subset of elements, which are to be
%                  computed are given in ind_local and the values produced are
%                  'INC_local'.  A time-step can then be realized by
%                  @code
%                  NU_local = U_local_ext(ind_local) - dt * INC_local
%                  @endcode example: fv_conv_explicit_space()
%
% optional fields of model:
% data_const_in_time : if this optional field is 1, the time evolution is
%                 performed with constant operators, i.e. only the
%                 initial-time-operators are computed and used throughout the
%                 time simulation.

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


% Martin Drohmann 9.12.2007 based on detailed_nonlin_evol_simulation.m by
% Bernard Haasdonk

if nargin ~= 2
  error('wrong number of parameters!');
end;

if model.verbose >= 9
  disp('entered detailed simulation ');
end;

if isempty(model_data)
  model_data = gen_model_data(model);
end;

% for a detailed simulation we do not need the affine parameter
% decomposition
model.decomp_mode = 0;

% initial values by midpoint evaluation
U0 = model.init_values_algorithm(model,model_data);

U = zeros(length(U0(:)),model.nt+1);
U(:,1) = U0(:);

if ~isfield(model,'data_const_in_time')
  model.data_const_in_time = 0;
end

if ~isfield(model, 'newton_regularisation')
  model.newton_regularisation = 0;
end

fl = cell(model.nt, 1);
% loop over fv-steps
if model.fv_impl_diff_weight ~= 0 && ~model.newton_solver
  % diffusive part is computed implicitly.
  clear_all_caches;
  [L_diff_I, const_diff_I] = model.implicit_operators_algorithm(model, model_data, []);
  clear_all_caches;
end

model.dt = model.T / model.nt;

nmax = 0;
nlast = 1;
Unewton = U0;
for t = 1:model.nt
  model.t = (t-1)*model.dt;
  model.tstep = t;
  if model.verbose > 19
    disp(['entered time-loop step ',num2str(t), ' for non linear evolution equation']);
  elseif model.verbose > 9
    fprintf('.');
  end;

  [inc, fluxes] = model.L_E_local_ptr(model, model_data, U(:,t), []);

  if model.debug == 1
    fl{t} = fluxes;
  end
  if all([model.fv_impl_conv_weight, model.fv_impl_diff_weight, ...
          model.fv_impl_react_weight] == [0, 0, 0])
    % purely explicit scheme
    U(:,t+1) = U(:,t) - model.dt * inc;
  elseif model.newton_solver
    V = zeros(length(U0(:)),model.newton_steps+1);
    exprhs = model.dt * inc;
    V(:,1) = U(:,t) - exprhs;
    n=0;
    residual=inf;
%    while true
%    n=n+1;
    for n=1:model.newton_steps;
      if model.verbose > 10
        disp(['Computing time step no. ', num2str(t), ...
              ', Newton step no. ', num2str(n), ...
              '... Residual norm: ']);
      end
      model.t = model.t + model.dt;
      Urhs     = model.L_I_local_ptr(model, model_data, V(:,n), []);
      [gradL_I, gL_I_offset] = ...
        model.implicit_gradient_operators_algorithm(model, model_data, ...
                                                    V(:,n),[]);
      model.t = model.t - model.dt;
      gradL_I  = gradL_I * model.dt + speye(size(gradL_I));
      rhs      = - V(:,n) + U(:,t) - exprhs - model.dt * Urhs;
      VU       = gradL_I \ rhs;
      residual = sqrt(model.inner_product(model, model_data, VU, VU));
      Vnew  = V(:,n) + VU;
      if model.newton_regularisation
        [i] = find(Vnew < 0);
        if(~isempty(i))
%          epsilon = 0.99*min(-V(i,n)./VU(i));
%          %        epsilon = 1;
          disp(' == regularisation == ');
          Vnew(i) = 0.000;
          %Vnew = V(:,n) + epsilon * VU;
        end
      end
      V(:,n+1) = Vnew;
      nres = Vnew - U(:,t) + exprhs + model.dt * Urhs;
      nres = sqrt(model.inner_product(model, model_data, nres, nres));
      if model.verbose > 10
        disp([num2str(residual), ' newton residual: ', num2str(nres)]);
      end
      if nres < model.newton_epsilon
        nmax = max(nmax, n);
        break;
      end
    end
    U(:,t+1) = V(:,n+1);
    Unewton = [Unewton, V(:,1:n+1)];
    nlast = nlast + n + 2;
  else % scheme with implicit parts, but without newton_scheme
    if ~model.data_const_in_time
      clear_all_caches;
      model.t=model.t+model.dt;
      [L_diff_I, const_diff_I] = model.implicit_operators_algorithm(model, model_data, []);
      model.t=model.t-model.dt;
      clear_all_caches;
    end
    inc = inc - const_diff_I;
    L_I = L_diff_I;
    L_I = L_I * model.dt + speye(size(L_I));
%    L_I = speye(size(L_diff_I));
%    keyboard;

    rhs = U(:,t) - model.dt * inc; % + model.dt * L_diff_I * U(:,t);
%    U(:,t+1) = rhs;
    U(:,t+1) = L_I \ rhs;
%    [U(:,t+1), flag] = bicgstab(L_I, rhs, 10e-5, 1000);
    if flag>0
      disp(['error in system solution, solver return flag = ', ...
        num2str(flag)]);
      keyboard;
    end;
  end
end

if model.newton_solver && model.verbose > 0
  disp(['computation ready: We needed at most ', num2str(nmax), ' Newton steps']);
end

sim_data.U = U;

if model.newton_solver
  sim_data.Unewton = Unewton;
end

if model.verbose > 9
  fprintf('\n');
end;

