function simulation_data = lin_evol_rb_simulation(model, reduced_data)
%function simulation_data = lin_evol_rb_simulation(model, reduced_data)
%
% function, which performs a reduced basis online simulation for
% the parameter vector mu, which is assumed to be set in the model
% class, i.e. a previous model = model.set_mu(model,[...]) is assumed
% to have taken place.
%
% allowed dependency of data: Nmax, N, M, mu
% not allowed dependency of data: H
% allowed dependency of computation: Nmax, N, M, mu
% not allowed dependency of computation: H
% Unknown at this stage: ---
%
% Required fields of model as required by the numerical operators and
%  mu_names       : the cell array of names of mu-components and
%                   for each of these stringt, a corresponding
%                   field in model is expected.
%  T              : end-time of simulation
%  nt             : number of timesteps to compute
%  error_norm     : 'l2' or 'energy'
%  L_I_inv_norm_bound: upper bound on implicit operator inverse norm
%                  constant in time, a single scalar
%  L_E_norm_bound: upper bounds on explicit operator
%                  constant in time, a single scalar
%  energy_norm_gamma: gamma >= 0 defining the weight in the energy norm
%  coercivity_bound_name : name of the function, which can
%                  computes a lower bound on the coercivity, e.g.
%                  fv_coercivity_bound
%
% plus fields required by coercivity_bound_name.m
%
% optional fields of model:
%  data_const_in_time : if this flag is set, then only operators
%                       for first time instance are computed
%  name_output_functional: if this field is existent, then an
%                            output estimation is performed and error.estimtations
%   starting_tim_step:  in t-partition case this is the starting time step
%                       for the actual t-partitioni time step
%   stopping_time_step: in t-partition this is the stopping time step for
%                       the actual t-partition
%
%   optional fileds of reduced_data:
%       Delta0: initial error
%
% generated fields of simulation_data:
%   a: time sequence of solution coefficients, columns are
%      'a(:,k)' = `a^(k-1)`
%   Delta:time sequence of `L^2`-posteriori error estimates
%            'Delta(k)' = `\Delta^{k-1}_N`
%         or energy-norm-posterior error estimates
%            'Delta_energy(k)' = `\bar \Delta^{k-1}_N`
%         depending on the field "error_norm" in model.
%
% if model.name_output_functional is set, then additionally, a
% sequence of output estimates 's(U(:,))' and error bound Delta_s is returned
%
% see the ***_gen_reduced_data() routine for specifications of the
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
% Markus Dihlmann Feb 2011

if model.verbose >= 10
  disp(['performing simulation for mu = ', mat2str(get_mu(model))]);
end;

if (~isequal(model.error_norm,'l2')) && ...
      (~isequal(model.error_norm,'energy'))
  error('error_norm unknown!!');
end;

%preparation for t-partition
if ~isfield(model,'transition_model')
    model.transition_model = 0;
end

if isfield(model,'starting_time_step')
    t_ind_start = model.starting_time_step;
else
    t_ind_start = 0;
end

if isfield(model, 'stopping_time_step')
    t_ind_stop = model.stopping_time_step;
else
    t_ind_stop = model.nt;
end


%cgrid = grid_geometry(params);
if ~model.transition_model
    nRB = length(reduced_data.a0{1});
else
    nRB = size(reduced_data.LL_E{1},2);
end
a = zeros(nRB,t_ind_stop-t_ind_start+1);   % matrix of RB-coefficients
Delta = zeros(1,t_ind_stop-t_ind_start+1); % vector of error estimators
model.dt = model.T/model.nt;


% initial data projection and linear combination
model.decomp_mode = 2;
sa0 = rb_init_values(model,[]);
a0 = lincomb_sequence(reduced_data.a0,sa0);
a(:,1) = a0;

%initial error
if isfield(reduced_data,'Delta0')
    Delta(1) = reduced_data.Delta0;
    %disp(['initial error is ',num2str(Delta(1))]);
end

% loop over time steps: computation of a(t+1)==a^t
LL_I_is_speye = -1; % unknown state here
for t = (t_ind_start+1):t_ind_stop
  model.t = (t-1)*model.dt;
  if model.verbose >= 20
    disp(['entered time-loop step ',num2str(t)]);
  end;

  % linear combination of all quantities
  % only once, if data is const in time

  if (t==t_ind_start+1) || (~model.data_const_in_time)
    [sLL_I, sLL_E, sbb, sK_II, sK_IE, sK_EE, sm_I, sm_E, sm] = ...
        rb_operators(model,[]);
    LL_I = lincomb_sequence(reduced_data.LL_I,sLL_I);
    LL_E = lincomb_sequence(reduced_data.LL_E,sLL_E);
    bb = lincomb_sequence(reduced_data.bb,sbb);
    K_II = lincomb_sequence(reduced_data.K_II,sK_II);
    K_IE = lincomb_sequence(reduced_data.K_IE,sK_IE);
    K_EE = lincomb_sequence(reduced_data.K_EE,sK_EE);
    m_I = lincomb_sequence(reduced_data.m_I,sm_I);
    m_E = lincomb_sequence(reduced_data.m_E,sm_E);
    m = lincomb_sequence(reduced_data.m,sm);
    if model.data_const_in_time
      % inverse of implicit matrix! L_I' = Id, L_E' = L_I^(-1) *
      % L_E, bb = L_I^(-1) bb
%      ILL_I = inv(LL_I);
      LL_E = LL_I \ LL_E;
      bb = LL_I \ bb;
      LL_I = speye(size(LL_I));
      LL_I_is_speye = 1;
    elseif max(max(LL_I- speye(size(LL_I))))<1e-12
      LL_I_is_speye = 1;
    end;
  end;

  % solve LL_I* a(:,t+1) = LL_E * a(:,t) + bb
  rhs = LL_E * a(:,t-t_ind_start) + bb;
%  warning('this does not work for pure space operators (independent of delta t). In that case use the commented line below.');
  %rhs = a(:,t) + model.dt * (LL_E * a(:,t) + bb);
  % check whether pure explicit problem:
  if (LL_I_is_speye == 1)
    a(:,t+1-t_ind_start) = rhs;
  elseif (LL_I_is_speye == -1) && isequal(LL_I, speye(size(LL_I)))
    a(:,t+1-t_ind_start) = rhs;
  else % solve linear system
       %    disp('check symmetry and choose solver accordingly!');
       %    keyboard;
       % nonsymmetric solvers:
       %  [U(:,t+1), flag] = bicgstab(Li,rhs,[],1000);
       %  [U(:,t+1), flag] = cgs(Li,rhs,[],1000);
       % symmetric solver, non pd:
       % see bug_symmlq for a very strange bug: cannot solve identity system!
       % [U(:,t+1), flag] = symmlq(Li,rhs,[],1000);
       %  [U(:,t+1), flag] = minres(Li,rhs,[],1000);
       % symmetric solver, pd:
       % [a(:,t+1), flag] = pcg(LL_I,rhs,[],1000);
       [a(:,t+1-t_ind_start), flag] = bicgstab(LL_I,rhs,[],1000);
       if flag>0
         disp(['error in system solution, solver return flag = ', ...
               num2str(flag)]);
         keyboard;
       end;
  end;

  % compute error estimate recursively
  % Delta^k := |L_I^(k-1)|^-1 * (res_norm(k-1)+
  % |L_E^(k-1)|*Delta^(k-1))
  %    keyboard;

  res_norm_sqr =  ...
      a(:,t+1-t_ind_start)' * K_II * a(:,t+1-t_ind_start) ...
      -2* a(:,t+1-t_ind_start)' * K_IE * a(:,t-t_ind_start) ...
      + a(:,t-t_ind_start)' * K_EE * a(:,t-t_ind_start) ...
      + m ...
      - 2* a(:,t+1-t_ind_start)' * m_I ...
      + 2* a(:,t-t_ind_start)' * m_E;
  % ensure real estimators (could be complex due to numerics)
  if res_norm_sqr>=0
    res_norm = sqrt(res_norm_sqr);
  else
    res_norm = 0;
  end;

  if isequal(model.error_norm,'l2')
    Delta(t+1-t_ind_start) = model.L_I_inv_norm_bound * ...
        (res_norm + ...
         model.L_E_norm_bound * Delta(t-t_ind_start));
  else
    % only sum up the squares of the residuals, multiplication with coeff
    % and squareroots are taken at the end.
    Rkplus1sqr = max(res_norm_sqr/(params.dt.^2), 0);
    Delta(t+1-t_ind_start) = Delta(t-t_ind_start)+Rkplus1sqr;
  end;
end;

if isequal(model.error_norm,'energy')
  % now perform scaling and squareroot of energy-error-estimate:
  alpha = model.coercivity_bound_ptr(model);
  if max(model.L_E_norm_bound>1) || (alpha == 0)
    % no reasonable estimation possible
    Delta(:) = nan;
  else
    C = (sqrt(1-model.L_E_norm_bound.^2)+1) * 0.5;
    Coeff = model.dt /(4*alpha*C*(1-model.energy_norm_gamma * C));
    Delta = sqrt(Delta * Coeff);
  end;
end;

% return simulation result
simulation_data       = [];
simulation_data.a     = a;
simulation_data.Delta = Delta;
simulation_data.LL_I  = LL_I;
simulation_data.LL_E  = LL_E;

if isfield(model,'name_output_functional')
  if isequal(model.error_norm,'l2')
    s = a' * reduced_data.s_RB;
    Delta_s = reduced_data.s_l2norm * Delta;
    simulation_data.s = s(:);
    simulation_data.Delta_s = Delta_s(:);
  else
    error(['output estimation not implemented concerning energy-norm.' ...
           ' choose l2 as error_norm!']);
  end;
end;

%| \docupdate 
