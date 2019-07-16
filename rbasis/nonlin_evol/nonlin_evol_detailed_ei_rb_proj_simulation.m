function ei_rb_proj_data = nonlin_evol_detailed_ei_rb_proj_simulation(model, detailed_data)
%function ei_rb_proj_data = nonlin_evol_detailed_ei_rb_proj_simulation(model, detailed_data)
% computes the detailed numerical scheme as described by 'model' for the
% parameter set in the same data structure 'model' but using empirical
% interpolated operators and projecting on the reduced basis in each time step.
%
% This function performs a general time evolution and computs the solution
% DOF-vector 'ei_rb_proj_data.U(:,t)' for all times 't=1,...,nt+1' instead of
% the real explicit operator, the empirical interpolation operator given in
% detailed_data is used. After each timestep, a projection on the RB-space is
% performed.  The result of this function should be identical to a pure reduced
% simulation with reconstruction.
%
%
% Generated fields of ei_sim_data:
%  U:        Dof Vector of result for each time instance (matrix of size
%            'H x (nt+1)')
%  c:        This is a matrix of size 'M x (nt+1)' holding the coefficient
%            vectors of the empirical interpolation of the space operator
%            for each time instance.
%  a:        This is a matrix of size 'N x (nt+1)' holding the coefficient
%            vectors of the reduced basis projection results.
%
% required fields of model:
% T             : final time
% nt            : number of time-intervals until 'T', i.e. 'nt+1' solution
%                 slices are computed
% init_values_algorithm: name of function for computing the initvalues-DOF with
%                 arguments (grid) example: init_values_cog()
% implicit_operators_algorithm: name of function for computing the 'L_I',
%                 'b_I'-operators with arguments '(grid)' example: 
%                 fv_operators_implicit()
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
% M :              number of colateral basis vectors, that will be used for
%                  empirical inteprolation of the explicit operator
% N :              number of reduced basis vectors, that will be used for
%                  projection of the detailed simulations
% inner_product_matrix_algorithm: name of method that computes the matrix for
%                  inner product computation.
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


% Bernard Haasdonk 24.6.2006

if nargin ~= 2
  error('wrong number of parameters!');
end;

if model.verbose >= 9
  disp(['entered detailed simulation with interpolation ',...
        'of the explicit operator and rb-projection']);
end;

% prepare online-data for online-ei-components (e.g. varying M)
% offline_data = rb_offline_prep(detailed_data,model);
% M is assumed to be set
% reduced_data = rb_online_prep(offline_data,model);
for oi = 1:length(detailed_data.BM)
  M = min( model.M, size(detailed_data.BM{oi},2) );

  % prepare online-data for online-ei-components (e.g. varying M)
  % M is assumed to be set
  QM{oi}    = detailed_data.QM{oi}(:,1:M);
  BM{oi}    = detailed_data.BM{oi}(1:M,1:M);
  TM{oi}    = detailed_data.TM{oi}(1:M);
end

if ~isfield(model, 'enable_error_estimator')
  model.enable_error_estimator = 0;
end

% test of the effect of magic-point exchange
%TM = TM(randperm(model.M));

model.decomp_mode = 0;

A = detailed_data.W;
RB = detailed_data.RB(:,1:model.N);

% initial values by midpoint evaluation
U0 = model.init_values_algorithm(model, detailed_data);

U = zeros(length(U0(:)),model.nt+1);
a = zeros(size(RB,2),model.nt+1);
if model.enable_error_estimator
  reslog = zeros(1, model.nt);
end

a(:,1) = RB' * A * U0(:);
U(:,1) = RB * a(:,1);
model.dt = model.T/model.nt;
operators_required = 1;
if isempty(model.data_const_in_time)
  model.data_const_in_time = 0;
end

impl_CRB_ind = detailed_data.implicit_crb_index;
expl_CRB_ind = detailed_data.explicit_crb_index;

% loop over fv-steps
for t = 1:model.nt  
  if model.verbose > 19
    disp(['entered time-loop step ',num2str(t)]);
  elseif model.verbose > 9
    fprintf('.');
  end;
  % get matrices and bias-vector 
  if operators_required
    model.t = (t-1)*model.dt;
    model.tstep = t;
    if all([model.fv_impl_conv_weight, model.fv_impl_diff_weight, model.fv_impl_react_weight] == [0, 0, 0])
      % purely explicit scheme
      tmpmodel = model;
      tmpmodel.decomp_mode = 0;
      [L_I_space, b_I] = model.implicit_operators_algorithm(model, ...
                                                            detailed_data);
    elseif ~model.newton_solver % implicit scheme without newton steps
      % get implicit contributions
      clear_all_caches;
      [L_I_red, b_I_red] = ...
        model.implicit_operators_algorithm(tmpmodel, detailed_data, []);
      clear_all_caches;
      L_I_space = QM{impl_CRB_ind} * (BM{impl_CRB_ind} \ ...
                                      L_I_red(TM{impl_CRB_ind},:));

      b_I = QM{impl_CRB_ind} * (BM{impl_CRB_ind} \ b_I_red(TM{impl_CRB_ind}));
      L_I_space = sparse(L_I_space);
    else
      b_I = 0;
    end
    % project b_I
    a_I = RB' * A * b_I;

    if model.data_const_in_time
      operators_required = 0;
    end;
  end;
  model.tstep = t;

  % compute real explicit operator result, but only use the magic-point entries
  % and perform an empirical interpolation instead
  inc = model.L_E_local_ptr(model, detailed_data, U(:,t), []);

  inc_local = inc(TM{expl_CRB_ind});
  %ei_coefficients = BM \ inc_local;
  %ei_inc = detailed_data.QM(:,1:model.M)*ei_coefficients;
  ei_inc = QM{expl_CRB_ind} * (BM{expl_CRB_ind} \ inc_local);

  % projection on RB space
  a_inc = RB' * A * ei_inc;

  if model.newton_solver
    updelta   = zeros(length(U0), model.newton_steps+1);
    updelta_a = zeros(model.N, model.newton_steps+1);
    exprhs_a = model.dt * a_inc;
    updelta_a(:,1) = a(:,t) - exprhs_a;
    n = 0;
    nmax = -1;
    Unewton = [];
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
      Urhs = model.L_I_local_ptr(model, detailed_data, updelta(:,n), []);

      % project result of L_I onto CRB
      Urhs_local     = Urhs(TM{impl_CRB_ind});
      Urhs_ei_coeffs = BM{impl_CRB_ind} \ Urhs_local;
      Urhs_ei        = QM{impl_CRB_ind}(:,1:model.M) * Urhs_ei_coeffs;

      % project Urhs_ei onto RB
      Urhs_a         = RB' * A * Urhs_ei;

      % get derivative DL_I (interpolation is done only implicitly)
      [gradL_I_local, gL_I_offset_local] = ...
        model.implicit_gradient_operators_algorithm(model, detailed_data, ...
                                                    updelta(:,n),[]);

      model.t = model.t - model.dt;
%      gL_I_offset_ei = QM{impl_CRB_ind} * (BM{impl_CRB_ind} \ ...
%                       gL_I_offset_local(reduced_data.TM_local));

%      gradL_I_ei = QM{impl_CRB_ind} * (BM{impl_CRB_ind} \ gradL_I_local);
      gradL_I_ei = gradL_I_local;

      % project onto RB space
      gradL_I_rb = RB' * A * gradL_I_ei * RB;

      gradL_I_rb  = gradL_I_rb * model.dt + speye(size(gradL_I_rb));

      rhs_a = -updelta_a(:,n) + a(:,t) - exprhs_a - model.dt * Urhs_a;

      delta_a = gradL_I_rb \ rhs_a;

      delta = RB * delta_a;

      residual = sqrt(model.inner_product(model, detailed_data, delta, delta));
      if model.verbose > 10
        disp(num2str(residual));
      end
      updelta(:,n+1) = updelta(:,n) + delta;
      updelta_a(:,n+1) = updelta_a(:,n) + delta_a;
      if residual < model.newton_epsilon
        nmax = max(nmax, n);
        break;
      end
    end
    a(:,t+1) = updelta_a(:,n+1);
    U(:,t+1) = updelta(:,n+1);
    Unewton = [ Unewton, updelta(:,1:n+1) ];
  else
    L_I = L_I_space * model.dt + speye(size(L_I_space));
    %  ei_inc = RB * a_inc;

    a_rhs = a(:,t) - model.dt * a_inc - model.dt * a_I;

    %  disp('check rhs and inc!!');
    %  keyboard;

    if isequal(L_I, speye(size(L_I)))
      a(:,t+1) = a_rhs;
      %    U(:,t+1) = rhs;
      U(:,t+1) = RB * a(:,t+1);
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
      a(:,t+1) = RB' * A * L_I * RB \ a_rhs;
      U(:,t+1) = RB * a(:,t+1);
      %    [U(:,t+1), flag] = bicgstab(L_I,rhs,[],1000); 
      if flag>0
        disp(['error in system solution, solver return flag = ', ...
          num2str(flag)]);
        keyboard;
      end;
    end;
    %  plot_element_data(grid,U(:,t),params);

    %U = RB * a;

    if model.enable_error_estimator
      % compute `\delta t R^k` explicitly and store in reslog(t)
      residuum  = U(:,t) - U(:,t+1) - model.dt * ei_inc;
      reslog(t) = sqrt(residuum' * A * residuum);
      %    res2      = -model.dt * (ei_inc - RB * a_inc);
      %    norm(res2) - norm(residuum) < 10*eps
    end
  end

end;

if model.newton_solver && model.verbose > 0
  disp(['computation ready: We needed at most ', num2str(nmax), ' Newton steps']);
end

if model.verbose > 9
  fprintf('\n');
end;


ei_rb_proj_data.U = U;
ei_rb_proj_data.a = a;

if model.newton_solver
  ei_rb_proj_data.Unewton = Unewton;
end

if model.enable_error_estimator
  ei_rb_proj_data.reslog = reslog;
end

