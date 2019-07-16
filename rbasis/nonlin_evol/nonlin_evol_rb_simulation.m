function simulation_data = nonlin_evol_rb_simulation(model, reduced_data)
%function simulation_data = nonlin_evol_rb_simulation(model, reduced_data)
%
% function, which performs a reduced basis online simulation for
% the parameter vector mu, which is assumed to be set in the model 
% struct, i.e. a previous model = model.set_mu(model, [...]) is assumed
% to have taken place.
%
% allowed dependency of data: Nmax, Mmax, N, M, mu
% not allowed dependency of data: H
% allowed dependency of computation: Nmax, Mmax, N, M, mu
% not allowed dependency of computation: H
% Unknown at this stage: ---
% 
% Required fields of model as required by the numerical operators and
%  mu_names       : the cell array of names of mu-components and
%                   for each of these stringt, a corresponding
%                   field in model is expected.
%  T              : end-time of simulation
%  nt             : number of timesteps to compute
%  L_E_local_name : name of the local space-discretization operator 
%                   evaluation with syntax
%
%                   INC_local = L_local(U_local_ext, ind_local,        
%                                 grid_local_ext, model)       
%                   where *_local indicates results on a set of
%                   elements, *_local_ext includes information including
%                   their neighbours, i.e. the grid must be known
%                   also for the neighbours, the values of the
%                   previous timstep must be known on the
%                   neighbours. The subset of elements, which are
%                   to be computed are given in ind_local and the
%                   values produced are INC_local. 
%                   Subsequent time-step operator would be defined
%                   as   NU_local = U(ind_local) - dt * INC_local
%  rb_operators   : function pointer to a method assembling the decomposed
%                   discretization operators
%
% optional fields of model:
%  data_const_in_time : if this flag is set, then only operators
%                       for first time instance are computed
%
% generated fields of simulation_data:
%  a(:,k)     : coefficient vector of the reduced
%               simulation at time t^{k-1}
%
% see the rb_***_online_prep() routine for specifications of the
% fields of reduced_data

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


% Bernard Haasdonk 15.5.2007

if nargin ~= 2
  error('wrong number of arguments: should be 2');
end

%% compute a posteriori error estimator (c.f. diploma thesis Martin Drohmann
%   p.  26ff.)
if ~isfield(model,'enable_error_estimator')
  model.enable_error_estimator = 0;
  model.Mstrich = 0;
end

if model.verbose >= 10
  disp(['performing simulation for mu = ',mat2str(get_mu(model))]);
end;

if ~isfield(model, 'Mstrich')
  model.Mstrich = 0;
end

if model.enable_error_estimator && model.Mstrich == 0;
  error('If error estimator is enabled, model.Mstrich must be > 0');
end

% possibly shrink reduced_data detected by this routine
reduced_data = nonlin_evol_reduced_data_subset(model,reduced_data);

% set init data
model.decomp_mode = 2;

sa0 = model.rb_init_values(model,[]);
a0  = lincomb_sequence(reduced_data.a0,sa0);

% perform simulation loop in time
model.dt = model.T/model.nt;
a        = zeros(model.N,model.nt+1);   % matrix of RB-coefficients
a(:,1)   = a0;

% in case of velocity-file access, the correct filename must be set here
% the following is only relevant in case of use of a 
% model.use_velocitymatrix_file and filecaching mode 2
if model.filecache_velocity_matrixfile_extract == 2;
  % change global velocity file to MMax velocity file generated in offline
  model.velocity_matrixfile = ...
      ['M',num2str(model.M(1)),'_','gridpart_Mmax',num2str(model.Mmax),'_',...
       model.velocity_matrixfile];
end

%disp('temporary true mass matrix is considered');
%M = reduced_data.M;

reduced_data.CE = {};

crbInd_imp = reduced_data.implicit_crb_index;
crbInd_exp = reduced_data.explicit_crb_index;

Mstrich = model.Mstrich;
M       = model.M;
% M plus Mstrich
MM      = M + Mstrich;

if model.enable_error_estimator
  % Lipschitz constants
  if ~isfield(model, 'C_E')
    C_E = 1;
  else
    C_E = model.C_E;
  end
  if ~isfield(model, 'C_I')
    C_I = 1;
  else
    C_I = model.C_I;
  end
  % initialize error estimator
  Delta = 0;
  % initialize error estimator log (strictly increasing)
  Deltalog = zeros(1, model.nt-1);
  % initialize log of residuum for each time step
  reslog   = zeros(1, model.nt-1);
  % initialize part for ei estimation for each time step
  eilog    = zeros(1, model.nt-1);
  % initialize part for newton error for each time step
  newtonlog    = zeros(1, model.nt);
end

if model.separate_CRBs
%  Mimpl = min([M+Mstrich, length(reduced_data.TM_local{crbInd_imp})]);
  Mexpl = cellfun(@length, reduced_data.TM_local(crbInd_exp,:));
  Mexpl = arrayfun(@(a,b) min(a,b), MM, Mexpl);
else
  Mimpl = cellfun(@length, reduced_data.TM_local(crbInd_exp,:));
  Mimpl = arrayfun(@(a,b) min(a,b), MM, Mimpl);
  Mexpl = Mimpl;
end

% model for explicit part (has different affine decomposition)
emodel             = model;
emodel.decomp_mode = 0;
emodel.M           = Mexpl;

if strcmp(model.RB_error_indicator, 'estimator') || ...
    strcmp(model.RB_error_indicator, 'ei_estimator_test')
  Mred_expl = min(M, Mexpl);
  Mred_impl = min(M, Mimpl);
else
  Mred_expl = Mexpl;
  Mred_impl = Mimpl;
end

reduced_data.CE     = cell(size(reduced_data.DE));
reduced_data.CE_red = cell(size(reduced_data.DE));
for i=1:length(reduced_data.DE(:))
  reduced_data.CE{i} = reduced_data.DE{i} / ...
    reduced_data.BM{i};
  reduced_data.CE_red{i} = reduced_data.DE{i}(:,1:Mred_impl(i)) ...
      / reduced_data.BM{i}(1:Mred_impl(i), 1:Mred_impl(i));
end

if ~isfield(reduced_data, 'time_split_map')
  reduced_data.time_split_map = [(1:model.nt+1)', ones(model.nt+1,1)];
end

%ll_E = zeros(Mexpl, model.nt);   % matrix of RB-coefficients
% loop over time steps: computation of a(t+1)==a^t
for t = 1:model.nt % either 1 or model.nt
  model.t      = (t-1)*model.dt;
  emodel.t     = model.t;
  emodel.tstep = t;
  if model.verbose >= 10
    disp(['entered time-loop step ',num2str(t)]);
  end;

  ti = reduced_data.time_split_map(t,2);
  tii = reduced_data.time_split_map(t+1,2);

  % asssemble affine parameter compositions of implicit parts
  % only once, if data is const in time
  if ~model.newton_solver && ((t==1) || (~model.data_const_in_time))
    if model.implicit_nonlinear
      reduced_data.grid = reduced_data.grid_local_ext{crbInd_imp,tii};
      tmpmodel = model;
      tmpmodel.decomp_mode = 0;
      clear_all_caches;
      [L_I_red, b_I_red] = ...
        model.implicit_operators_algorithm(tmpmodel, reduced_data,...
                                           reduced_data.TM_local{crbInd_imp,tii});
      clear_all_caches;

      LL_I_red = L_I_red * reduced_data.RB_local_ext{crbInd_imp,tii};
      LL_I     = reduced_data.CE{crbInd_imp,tii}(:,1:Mred_impl(tii)) * ...
                    LL_I_red;

      bb_I     = reduced_data.CE{crbInd_imp,tii}(:,1:Mred_impl(tii)) * ...
                   b_I_red;

    else
      [sLL_I, sbb_I] = ...
        model.rb_operators(model,[]); % detailed_data omitted
      if ~isempty(reduced_data.LL_I)
        LL_I = lincomb_sequence(reduced_data.LL_I,sLL_I);
      else
        LL_I = sparse([],[],[],model.N,model.N);
      end;
      if ~isempty(reduced_data.bb_I)
        bb_I = lincomb_sequence(reduced_data.bb_I,sbb_I);
      else
        bb_I = sparse([],[],[],model.N,1);
      end;
    end;
  end;

  % setup linear system     A * a(:,t+1) = rhs
  % with A = Id + dt * LL_I
  % rhs = a(:,t) - dt * BB_I - dt * CE * ll_E(a(:,t))
  %                  where ll_E is the local operator evaluation

  if ~model.newton_solver
    A = LL_I * model.dt + speye(size(LL_I));
  else
    A = 0;
    bb_I = 0;
  end

  % determine local values by reconstructing U values from reduced
  % basis parts
  U_local_ext = reduced_data.RB_local_ext{crbInd_exp,ti} * a(:,t);

  % the following performs
  %  U_local = L_local(U_local_ext, TM_local, [], ...
  %		    reduced_data.grid_local_ext);
  %
  % note the use of emodel
  reduced_data.grid = reduced_data.grid_local_ext{crbInd_exp,ti};

  le_local_inc = emodel.L_E_local_ptr(emodel, reduced_data, ...
                                      U_local_ext, ...
                                      reduced_data.TM_local{crbInd_exp,ti});
%  le_local_inc(U_local_ext(reduced_data.TM_local{crbInd_exp,tii}) + le_local_inc > 1) = 0;

  if model.newton_solver
    nmax = -1;
    updelta_a = zeros(model.N, model.newton_steps+1);
    updelta_a(:,1) = a(:,t);
    n=0;

    for n=1:model.newton_steps
      if model.verbose > 10
        disp(['Computing time step no. ', num2str(t), ...
              ', Newton step no. ', num2str(n), ...
              '... Residual norm: ']);
      end
      model.decomp_mode = 0;
      updelta_rec = reduced_data.RB_local_ext{crbInd_imp,tii} * updelta_a(:,n);
      if model.debug && any(updelta_rec < -10*eps)
        disp(num2str(min(updelta_rec)));
        warning('RBMatlab:rb_simulation:runtime', ...
                'reconstructed solutions at interpolation points is negative!');
      end
      model.t = model.t + model.dt;
      % updelta_rec(updelta_rec < 0) = 0.001;
      reduced_data.grid = reduced_data.grid_local_ext{crbInd_imp,tii};
      li_local_inc = model.L_I_local_ptr(model, reduced_data, ...
                                         updelta_rec, ...
                                         reduced_data.TM_local{crbInd_imp,tii});

      [gradL_I, gL_I_offset] = ...
        model.implicit_gradient_operators_algorithm(model, reduced_data, ...
                                          updelta_rec,...
                           reduced_data.TM_local{crbInd_imp,tii}(1:Mred_impl(tii)));
      model.t = model.t - model.dt;
      if any(gL_I_offset ~= 0)
        warning('offset has non-zero entries!');
      end

      gradL_I  = reduced_data.CE_red{crbInd_imp,tii} * gradL_I ...
                   * reduced_data.RB_local_ext{crbInd_imp,tii};

      gradL_I  = gradL_I * model.dt + speye(size(gradL_I));

      opts.LT = true;

      if ti ~= tii
        theta_update_expl = linsolve(reduced_data.BM{crbInd_exp,ti}(1:MM(ti),1:MM(ti)), ...
                                    le_local_inc(1:MM(ti)), ...
                                    opts);
        theta_update_impl = linsolve(reduced_data.BM{crbInd_imp,tii}(1:MM(tii),1:MM(tii)), ...
                                    li_local_inc(1:MM(tii)), ...
                                    opts);

%        minlength = min(length(theta_update_impl), length(theta_update_expl));
%        theta_update = theta_update_expl(1:minlength) + theta_update_impl(1:minlength);

        arhs = reduced_data.DE{crbInd_exp,ti}(:,1:Mred_expl(ti)) ...
                          * theta_update_expl(1:Mred_expl(ti)) ...
               + reduced_data.DE{crbInd_imp,tii}(:,1:Mred_impl(tii)) ...
                          * theta_update_impl(1:Mred_impl(tii));
      else
        theta_update = linsolve(reduced_data.BM{crbInd_exp,ti}(1:MM(ti),1:MM(ti)), ...
                                le_local_inc(1:MM(ti)) + li_local_inc(1:MM(ti)), ...
                                opts);

        arhs = reduced_data.DE{crbInd_exp,ti}(:,1:Mred_expl(ti)) ...
                          * theta_update(1:Mred_expl(ti));
      end

      full_rhs = - updelta_a(:,n) + a(:,t) - model.dt * arhs;

      delta_a = gradL_I \ full_rhs;

      updelta_a(:,n+1) = updelta_a(:,n) + delta_a;
      nres_a = updelta_a(:,n+1) -a(:,t) + model.dt * arhs;
      nres = sqrt(nres_a' * nres_a);
      if model.verbose > 10
        disp(num2str(nres));
      end
      if nres < model.newton_epsilon
        nmax = max(nmax, n);
        break;
      end
    end
    a(:,t+1) = updelta_a(:,n+1);
  else
    % compute coefficients of update in CRB Space `\W_{M+M'}`.
    opts.LT = true;
    theta_update_expl = linsolve(reduced_data.BM{crbInd_exp,ti}(1:MM(ti),1:MM(ti)), ...
                                 le_local_inc, opts);
    theta_update_impl = 0;

    rhs = a(:,t) - model.dt * bb_I ...
                 - model.dt * reduced_data.DE{crbInd_exp,ti}(:,1:Mred_expl(ti)) ...
                            * theta_update_expl(1:Mred_expl(ti));

    % check if an LES is to be solved, or problem is an explicit problem
    if isequal(A, speye(size(A)))
      a(:,t+1) = rhs;
    else
      % solve linear system without Newton solver
%        disp('check symmetry and choose solver accordingly!');
%       keyboard;
       % nonsymmetric solvers:
       [a(:,t+1), flag, relres, iter] = bicgstab(A,rhs,1e-10,1000);
       % [a(:,t+1), flag] = cgs(Li,rhs,[],1000); 
       % symmetric solver, non pd:
       % see bug_symmlq for a very strange bug: cannot solve identity system!
       % [a(:,t+1), flag] = symmlq(Li,rhs,[],1000);
       %  [a(:,t+1), flag] = minres(Li,rhs,[],1000); 
       % symmetric solver, pd:
       %[a(:,t+1), flag] = bicgstab(A,rhs,[],1000); 
       if(t==1 && model.debug)
         disp(cond(A));
         disp(all(A'==A));
         disp(eig(full(A)));
       end
       if(t==1) % && model.verbose > 10)
         disp(['relres=', num2str(relres), ' iter=', num2str(iter)]);
       end
       % exact solver:
       % a(:,t+1)=A \ rhs;
       if flag>0
         disp(['error in system solution, solver return flag = ', ...
               num2str(flag)]);
         keyboard;
       end
    end
  end

  if model.enable_error_estimator
    if tii > 2
      error('error estimator does not work with varying crbs for different time slices yet');
    end
    ainc       = a(:,t) - a(:,t+1);
    theta_update   = theta_update_expl + theta_update_impl;

    if model.newton_solver
      lipschitzConst = C_I^(model.nt - t + 2)*C_E^(model.nt - t + 1);
      newtonlog(t)   = nres + model.newton_epsilon;
    else
      lipschitzConst = C_E^(model.nt - t + 1);
    end
    % compute the residuum `\delta t norm(R^k)`
    % TODO revise the following 
    residuum  = sqrt( abs( ainc' * ainc ...
                      - 2*model.dt * ainc' ...
                                   * ( reduced_data.DE{crbInd_exp,ti}(:,1:M(ti)) ...
                                       * theta_update(1:M(ti)) ...
                                       + reduced_data.DE{crbInd_imp,ti}(:,1:M(ti)) ...
                                       * theta_update(1:M(ti)) ) ...
                      + model.dt^2 * theta_update(1:M(ti))' ...
                                   * reduced_data.Mmass{1,ti}(1:M(ti),1:M(ti)) ...
                                   * theta_update(1:M(ti))) );
    reslog(t) = residuum;

    % compute the empirical interpolation error estimation
    eilog(t) = model.dt*lipschitzConst ...
                     * sqrt(abs( ...
                     theta_update(M(tii)+1:MM(tii))' ...
                        * reduced_data.Mmass{1,tii}(M(tii)+1:MM(tii),M(tii)+1:MM(tii)) ...
                     * theta_update(M(tii)+1:MM(tii)) ));

    Delta    = Delta + eilog(t) ...
                     + lipschitzConst*newtonlog(t) ...
                     + lipschitzConst*residuum;
    if ~isreal(Delta)
      warning('RBMatlab:rb_simulation:runtime',...
              'error estimator returned complex result');
      Delta = real(Delta);
    end
    Deltalog(t) = Delta;
  end
end; % timeloop

if model.newton_solver && model.verbose > 0
  disp(['computation ready: We needed at most ', num2str(nmax), ' Newton steps']);
end

if model.enable_error_estimator
  simulation_data.Delta  = [0, Deltalog];
  simulation_data.reslog = [0, reslog];
  simulation_data.eilog  = [0, eilog];
end

% prepare return parameters
simulation_data.a = a;
%| \docupdate
