function ei_sim_data = nonlin_evol_detailed_ei_simulation(model, detailed_data)
%function ei_sim_data = nonlin_evol_detailed_ei_simulation(model, detailed_data)
% computes the detailed numerical scheme as described by 'model' for the
% parameter set in the same data structure 'model' but using empirical
% interpolated operators.
%
% This function performs a general time evolution and computs the solution
% DOF-vector 'ei_sim_data.U(:,t)' for all times 't=1,...,nt+1' instead of the
% real explicit operator, the empirical interpolation operator given in
% detailed_data is used.  So this function can be used for testing the
% stability of the empirical interpolation operator.
%
% Generated fields of ei_sim_data:
%  U:        Dof Vector of result for each time instance (matrix of size
%            'H x (nt+1)')
%  c:        This is a matrix of size 'M x (nt+1)' holding the coefficient
%            vectors of the  the empirical interpolation of the space operator
%            for each time instance.
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
        'of the explicit operator']);
end;

if ~iscell(detailed_data.BM)
  detailed_data.BM = { detailed_data.BM, detailed_data.BM, detailed_data.BM };
  detailed_data.QM = { detailed_data.QM, detailed_data.QM, detailed_data.QM };
  detailed_data.TM = { detailed_data.TM, detailed_data.TM, detailed_data.TM };
end

M = min(model.M, max(cellfun(@(X)(size(X,2)), detailed_data.BM)));

if model.separate_CRBs
  Mexpl = min(model.M, size(detailed_data.BM{detailed_data.explicit_crb_index},2));
else
  Mexpl = M;
end
QM = cell(1,length(detailed_data.BM));
BM = cell(1,length(detailed_data.BM));
TM = cell(1,length(detailed_data.BM));
BMinv = cell(1,length(detailed_data.BM));

for oi=1:length(detailed_data.BM)
  M = min(model.M, max(size(detailed_data.BM{oi},2)));
  % prepare online-data for online-ei-components (e.g. varying M)
  QM{oi}    = detailed_data.QM{oi}(:,1:M);
  BM{oi}    = detailed_data.BM{oi}(1:M,1:M);
  TM{oi}    = detailed_data.TM{oi}(1:M);
  BMinv{oi} = inv(BM{oi});
end

% test of the effect of magic-point exchange
%TM = TM(randperm(model.M));

model.decomp_mode = 0;

% initial values by midpoint evaluation
U0 = model.init_values_algorithm(model, detailed_data);

U = zeros(length(U0(:)),model.nt+1);
U(:,1) = U0(:);

c = zeros(Mexpl,model.nt);

model.dt = model.T/model.nt;
operators_required = 1;
if isempty(model.data_const_in_time')
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
      [LL_I, bb_I] = ...
        model.implicit_operators_algorithm(tmpmodel, detailed_data, []);
    elseif ~model.newton_solver % implicit scheme without newton steps
      tmpmodel = model;
      tmpmodel.decomp_mode = 0;
      clear_all_caches;
      [L_I_red, b_I_red] = ...
        model.implicit_operators_algorithm(tmpmodel, detailed_data, []);
      clear_all_caches;
      LL_I = QM{impl_CRB_ind} * BMinv{impl_CRB_ind} ...
               * L_I_red(TM{impl_CRB_ind},:);
      bb_I = QM{impl_CRB_ind} * BMinv{impl_CRB_ind} ...
               * b_I_red(TM{impl_CRB_ind});
      LL_I = sparse(LL_I);
    end;
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

  c(:,t) = BM{expl_CRB_ind} \ inc_local;

  ei_inc = QM{expl_CRB_ind} * (BM{expl_CRB_ind} \ inc_local);

  if model.newton_solver
%    dTM   = detailed_data.TM{impl_CRB_ind};
%    mask = zeros(1, detailed_data.grid.nelements);
%    nbi  = detailed_data.grid.NBI(dTM,:);
%    i    = find(nbi > 0);
%    % get indices of TM's neighbour-neighbours
%    nnbi            = detailed_data.grid.NBI(unique(nbi(i)),:);
%    ni              = find(nnbi > 0);
%    [nnbuniq,tally] = unique(sort(nnbi(ni)),'first');
%    tally           = [tally(2:end);length(nnbi(ni))+1] - tally;
%    mask(nnbuniq)   = tally;
%    mask(nbi(i))    = 5;
%    mask(dTM)        = 6;
%    mask(mask < 2)  = 0;
%
%    eind            = find(mask);
%    eind_speye = sparse(1:length(eind), eind, ones(1,length(eind)), ...
%                        length(eind), ...
%                        detailed_data.grid.nelements);
%
%    reduced_data.grid = gridpart(detailed_data.grid,eind);
%
%    glob2loc              = zeros(detailed_data.grid.nelements,1);
%    glob2loc(eind)        = 1:length(eind);
%    reduced_data.TM_local = glob2loc(dTM);

    updelta = zeros(length(U0),model.newton_steps+1);
    exprhs = model.dt * inc;
    updelta(:,1) = U(:,t) - exprhs;
    n = 0;
    nmax = -1;
    Unewton = U0;
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

      % get derivative DL_I (interpolation is done only implicitly)
      [gradL_I_local, gL_I_offset_local] = ...
        model.implicit_gradient_operators_algorithm(model, detailed_data, ...
                                                    updelta(:,n),[]);

      model.t = model.t - model.dt;
%      gL_I_offset_ei = QM{impl_CRB_ind} * (BM{impl_CRB_ind} \ ...
%                       gL_I_offset_local(reduced_data.TM_local));

%      gradL_I_ei = QM{impl_CRB_ind} * (BM{impl_CRB_ind} \ gradL_I_local);
      gradL_I_ei = gradL_I_local;

      gradL_I  = gradL_I_ei * model.dt + speye(size(gradL_I_ei));

      rhs = - updelta(:,n) + U(:,t) - exprhs - model.dt * Urhs_ei;

      delta = gradL_I \ rhs;

      % project delta on CRB
      delta_local     = delta(TM{impl_CRB_ind});
      delta_ei_coeffs = BM{impl_CRB_ind} \ delta_local;
      delta_ei        = QM{impl_CRB_ind}(:,1:model.M) * delta_ei_coeffs;

      residual = sqrt(model.inner_product(model, detailed_data, delta_ei, delta_ei));
      if model.verbose > 10
        disp(num2str(residual));
      end
      updelta(:,n+1) = updelta(:,n) + delta_ei;
      iltz = find(updelta(:,n+1) < 0);
      if ~isempty(iltz)
          updelta(iltz,n+1) = 0.001;
      end
      if residual < model.newton_epsilon
        nmax = max(nmax, n);
        break;
      end
    end
    U(:,t+1) = updelta(:,n+1);
    Unewton = [ Unewton, updelta(:,1:n+1) ];
  else
    L_I = LL_I * model.dt + speye(size(LL_I));

    rhs = U(:,t) - model.dt * ei_inc - model.dt * bb_I;

    %  disp('check rhs and inc!!');
    %  keyboard;

    if isequal(L_I, speye(size(L_I)))
      U(:,t+1) = rhs;
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
      U(:,t+1) = L_I \ rhs;
      %    [U(:,t+1), flag] = bicgstab(L_I,rhs,[],1000);
      if flag>0
        disp(['error in system solution, solver return flag = ', ...
              num2str(flag)]);
        keyboard;
      end;
    end;
  end
%  plot_element_data(grid,U(:,t));

end;

if model.newton_solver && model.verbose > 0
  disp(['computation ready: We needed at most ', num2str(nmax), ' Newton steps']);
end

if model.verbose > 9
  fprintf('\n');
end;

ei_sim_data.U = U;
ei_sim_data.c = c;
if model.newton_solver
  ei_sim_data.Unewton = Unewton;
end

