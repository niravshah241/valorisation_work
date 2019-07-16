function ei_sim_data = nonlin_evol_detailed_local_ei_simulation(model, detailed_data)
%function ei_sim_data = nonlin_evol_detailed_local_ei_simulation(model, detailed_data)
%
% function performing a general time evolution and computing the
% solution DOF-vector ei_sim_data.U(:,t) for all times t=1,...,nt+1
% instead of the real explicit operator, the empirical
% interpolation operator given in detailed_data is used.
% So this function can be used for testing the stability of the
% empirical interpolation operator. The coefficient vectors of the
% empirical interpolation of the space operator are returned as
% columns of ei_sim_data.c
% 
% In contrast to the nonlin_evol_detailed_ei_simulation, here, the
% ei_operator is evaluated locally.
%
% required fields of model:
% T             : final time
% nt            : number of time-intervals until T, i.e. nt+1
%                 solution slices are computed
% init_values_algorithm: name of function for computing the
%                 initvalues-DOF with arguments (grid)
%                 example: init_values_cog
% implicit_operators_algorithm: name of function for computing the
%                 L_I, b_I-operators with arguments (grid)
%                 example fv_operators_implicit
% data_const_in_time : if this optional field is 1, the time
%                 evolution is performed with constant operators,
%                 i.e. only the initial-time-operators are computed
%                 and used throughout the time simulation.
% L_E_local_ptr :  function pointer to the local space-discretization
%                  operator evaluation with syntax
%
%                   INC_local = L_local(U_local_ext, ind_local,
%                                      grid_local_ext)
%                   where *_local indicates results on a set of
%                   elements, *_local_ext includes information including
%                   their neighbours, i.e. the grid must be known
%                   also for the neighbours, the values of the
%                   previous timstep must be known on the
%                   neighbours. The subset of elements, which are
%                   to be computed are given in ind_local and the
%                   values produced are INC_local.
%                   A time-step can then be realized by
%                        NU_local = U_local_ext(ind_local) - dt * INC_local
% M : number of colateral basis vectors, that will be used for
%                   empirical inteprolation of the explicit operator

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

%grid = detailed_data.grid;

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
    if isfield(model,'model_type') && ...
	  isequal(model.model_type, 'implicit_nonaffine_linear')
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
    else
      tmpmodel = model;
      tmpmodel.decomp_mode = 0;
      [LL_I, bb_I] = ...
        model.implicit_operators_algorithm(tmpmodel, detailed_data, []);
    end;
    if model.data_const_in_time
      operators_required = 0;
    end;
  end;
  model.tstep = t;
  L_I = LL_I * model.dt + speye(size(LL_I));

  % compute real explicit operator result, but only use the
  % magic-point entries and perform an empirical interpolation instead
  %inc = model.L_E_local_ptr(model, detailed_data, U(:,t), []);
  %inc_local = inc(TM{expl_CRB_ind});

  % local evaluation:

  %  keyboard;
  
  grid = detailed_data.grid;
  eind = TM{1};
  eind_ext = index_ext(grid,...
		       eind,...
		       model.local_stencil_size);
  U_local_ext = U(eind_ext,t);
  %  exp_ind_ext =  ...
  %  exp_ind_local = ...
  %  U_local_ext = ...
  tmp = zeros(1,grid.nelements);
  tmp(eind_ext) = 1:length(eind_ext);
  eind_local = tmp(eind)
  
  if ~isfield(detailed_data,'grid_local_ext');
    detailed_data.grid_local_ext = {...
	gridpart(grid,eind_ext)};
  end;
  
  inc_local = model.L_E_local_ptr(model, detailed_data, ...
				  U_local_ext, eind_local);

  %  ei_coefficients = BM \ inc_local;
  %  ei_inc = detailed_data.QM(:,1:model.M)*ei_coefficients;
  %  c(:,t) = BM{expl_CRB_ind} \ inc_local;
  
  ei_inc = QM{expl_CRB_ind} * (BM{expl_CRB_ind} \ inc_local);

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
%  plot_element_data(grid,U(:,t));

end;
if model.verbose > 9
  fprintf('\n');
end;

ei_sim_data.U = U;
ei_sim_data.c = c;
%| \docupdate 
