function [LU_fnames] = ei_operator_collect_files(model, model_data, Mtrain, params)
%function [LU_fnames] = ei_operator_collect_files(model,model_data,Mtrain,params)
% collects operator evaluations on a sample of snapshots
%
% This function collects operator evaluations applied to sets of snapshots.
% Given a list of operators parameter dependent operators `L_q(\mu) \in {\cal
% W}_h, q=1,\ldots,Q` which are `H`-independent Dof dependent as described in
% [DOH10] and a set of parameters `M \subset {\cal P} \subset {\cal R}^p`, this
% function computes operator evaluations `L_q(\mu)[u_h(\mu;t^k)]` for all `\mu
% \in M` and `t^k=k\Delta t`, and stores the results in matlab files.
%
% parameters:
%      Mtrain     : a matrix specifying the parameter set `M` for which
%        detailed simulations shall be computed. The column vectors of the
%        matrix each contain a parameter vector.
%      params     : a structure controlling the behaviour of this function
%
% required fields of model:
%      collect_newton_steps : boolean flag, indicating whether the intermediate
%        newton steps computed by a detailed simulaton with a Newton scheme
%        shall also be considered when applying the operators 
%      ei_detailed_savepath : path, where the trajectories computed as results
%        of detailed simulations are stored. This files must contain a
%        structure named 'sim_data' with a field 'U' of dimesion
%        '\#DOFs x model.nt +1' or a field 'Unewton' of dimension
%        '\#DOFs x model.nt + 1 + \#intermediate Newton steps' in case of a
%        Newton scheme and 'model.collect_newton_steps==true'. The files are
%        created by save_detailed_simulations() for all parameters given by
%        'Mtrain'
%      ei_operator_savepath : path, where the LU-files are stored.  These files
%        each contain an array named 'LU' of the same dimension as 'sim_data.U'
%        respectively 'sim_data.Unewton' from the files stored in
%        'ei_detailed_savepath'.
%      separate_CRBs        : boolean flag, indicating whether, the operator
%        evaluations shall be stored in different directories for the later
%        computation of separate collateral reduced basis spaces for each
%        operator, or combined a single directory.
%
% optional fields of model:
%   force_delete : if this option is set to true (the default), directories
%                  with data inconsistent with the current settings are
%                  deleted.
%
% required fields of params:
%      ei_space_operators: cell array of function pointers to spacial operators
%        `L_q, q=1,\ldots,Q` for which evaluations on snapshots shall be
%        computed and stored.
%
% optional fields of params:
%      num_cpus:           number of cpu cores that can be used for parallel
%        computations. This parameter has only an effect, when Matlab's
%        parallel computation capabilities are started by initiating a
%        matlabpool.  (default == 4)

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


% paths are to be understood relative to $RBMATLABTEMP

% Bernard Haasdonk 15.8.2007

% for each mu in Mtrain compute exact detailed trajectory u_H^k(mu)
npar = size(Mtrain,2);

if (model.verbose > 0)
  disp(['generating detailed data in ',...
       fullfile(rbmatlabtemp,model.ei_detailed_savepath)]);
end;

save_detailed_simulations(model, model_data, Mtrain, model.ei_detailed_savepath);

skip_computation = 0;
savepath = fullfile(getenv('RBMATLABTEMP'),...
                    model.ei_operator_savepath);

if (model.verbose > 0)
  disp(['generating ei_operator snapshot data in ',...
        fullfile(rbmatlabtemp,model.ei_operator_savepath)]);
end;

if ~isfield(model, 'force_delete')
  model.force_delete = 1;
end

if ~exist(savepath,'dir');
  [s,m,id] =mkdir(getenv('RBMATLABTEMP'), ...
                          model.ei_operator_savepath);
  if ~s
    error('problem in creating directory!');
  end;
elseif exist(fullfile(savepath,'settings.mat'), 'file')
  % checking consistency of existing
  disp('directory exists, checking for correct data.');
  tmp = load(fullfile(savepath,'settings.mat'));

  if ~isequal(Mtrain, tmp.Mtrain)
    if ~model.debug && model.force_delete
      disp(['Mtrain in precomputed data and current are inconsistent!',...
            ' I am deleting the old data now!']);
      delete(fullfile(savepath,'*.mat'))
      ls(savepath);
    else
      error(['Mtrain in precomputed data and current',...
             ' simulation are inconsistent! Please delete and restart!']);
    end
  else
    % fields to skip in comparison
    ignorelist = {'Mmax','ei_Mmax','RB_stop_epsilon','RB_stop_Nmax',...
                 'ei_target_error', 'velocity_coefficients_ptr','Mstrich',...
                  model.mu_names{:},'force_delete','RB_error_indicator',...
                  'ei_time_splits','verbose', 'adaptive_time_split_mode', ...
                  'extend_crb', ...
                  'ei_stop_epsilon', 'ei_stop_epsilon_first', ...
                  'minimum_time_slice', 'time_split_Mmax'};
    ignorelist = [ignorelist, tmp.model.filecache_ignore_fields_in_model];
    [iseq, a,b,c] = structcmp(model,tmp.model,ignorelist);
    if ~iseq
      disp('fields of params and stored params differs:')
      disp('differing fields:');
      cellfun(@disp,a);
      disp('additional fields in params:');
      cellfun(@disp,b);
      disp('additional fields in stored params:');
      cellfun(@disp,c);
      if ~model.debug && model.force_delete
        disp(['Mtrain in precomputed data and current are inconsistent!',...
          ' I am deleting the old data now!']);
        delete(fullfile(savepath,'*.mat'))
        ls(savepath);
      else
        error(['parameters of precomputed data and current',...
               ' simulation are inconsistent! Please delete files at path ', ...
                savepath, ' and restart!']);
      end

    else
      if exist(fullfile(savepath,'UNFINISHED.lock'),'file');
        disp(['detected UNFINISHED.lock in target directory, ',...
              'resuming previous computations...']);
      else
        disp('skipping ei_operator evaluations on snapshots, as results existing')
        skip_computation = 1;
      end
    end
  end
else
  disp('stray directory found! I am cleaning this one, now!');
  delete(fullfile(savepath,'*.mat'))
  ls(savepath);
end

save(fullfile(savepath,'UNFINISHED.lock'),'savepath');
save(fullfile(savepath,'settings.mat'),'Mtrain','model');

LU_fnames = cell(length(params.ei_space_operators),npar);
%BU_fnames = cell(length(params.ei_space_operators),npar);

% for each of these trajectories,
% select time-instants and apply tobe-interpolated operator

model.dt = model.T / model.nt;

if ~iscell(params.ei_space_operators)
  params.ei_space_operators = { params.ei_space_operators };
end

if ~isfield(params, 'num_cpus')
  params.num_cpus = 4;
end
num_cpus = params.num_cpus;

%if ~isfield(params, 'separate_CRBs')
%  params.separate_CRBs = false;
%end
%

for oi = 1:length(params.ei_space_operators)
  for mout = 0:floor((npar-1)/num_cpus)
    LUc = cell(1, npar);
    parfor m = (mout*num_cpus+1):min(npar,(mout+1)*num_cpus)
%      mm = m - mout*num_cpus;
      tmp_params = params;
      opname = tmp_params.ei_space_operators{oi};
      tmp_model_data = model_data;
      tmp_model      = model;
      oiIdent = num2str(oi);
      disp(['processing parameter vector ',num2str(m),'/',num2str(npar)]);
      LU_fnames{oi,m} = fullfile(tmp_model.ei_operator_savepath,...
                                 ['LU_',oiIdent,'_',num2str(m),'.mat']);
%      BU_fnames{oi,m} = fullfile(tmp_model.ei_operator_savepath,...
%                                 ['BU_',oiIdent,'_',num2str(m),'.mat']);
      skip_file = exist(LU_fnames{oi,m},'file');% && exist(BU_fnames{oi,m},'file')
      if ~skip_computation && ~skip_file
        % determine initial size of vectors:
        sim_data = load_detailed_simulation(m,model.ei_detailed_savepath,model);
        if model.newton_solver && model.collect_newton_steps
          U = sim_data.Unewton;
        else
          U = sim_data.U;
        end
        TLU = zeros(size(U,1),0);
        %TBU = zeros(size(U,1),0);
        newmodel = model.set_mu(model, Mtrain(:,m));
        if isfield(model,'model_type') && isequal(newmodel.model_type, 'implicit_nonaffine_linear')
          clear_all_caches;
          [tmp_model_data.implicit_operator, tmp_model_data.implicit_constant] = ...
            model.operators_diff_implicit(newmodel, tmp_model_data, []);
          clear_all_caches;
        end
        if newmodel.newton_solver
          for tn = 1:size(U,2);
            fprintf('.');
            if newmodel.debug && ~newmodel.data_const_in_time
              error('time dependent data is not yet implemented for newton schemes');
            end
            [nu, bu] = opname(newmodel, tmp_model_data, U(:,tn), []);
            if ~isempty(bu) && tn==1
              TLU = [TLU, bu];
            end;
            TLU = [TLU, nu(:)];
          end
        else
          for tn = 1:length(newmodel.ei_time_indices);
            fprintf('.');
            ti = newmodel.ei_time_indices(tn);
            % note the time-shift by 1 in the following (ti=1 is initial data)
            newmodel.t = (ti-1)*newmodel.dt;
            newmodel.tstep = tn;
            [nu, bu] = opname(newmodel, tmp_model_data, U(:,ti), []);
            if ~isempty(bu) && tn==1
              %          TBU = [TBU, bu];
              TLU = [TLU, bu];
              %          LU  = TBU;
              %          save(fullfile(rbmatlabtemp,BU_fnames{oi,m}),'LU');
            end;
            TLU = [TLU, nu(:)];
          end
        end
        fprintf('\n');
        LUc{m} = TLU;
      end
    end
    if ~skip_computation
      for m = (mout*num_cpus+1):min(npar,(mout+1)*num_cpus)
        skip_file = exist(LU_fnames{oi,m},'file');% && exist(BU_fnames{oi,m},'file')
        if ~skip_file
          %        mm = m - mout*num_cpus;
          LU = LUc{m};
          save(fullfile(rbmatlabtemp,LU_fnames{oi,m}),'LU');
        end
      end
    end
  end
end

if ~model.separate_CRBs
  LU_fnames = reshape(LU_fnames, 1, size(LU_fnames,1)*size(LU_fnames,2));
end

if exist(fullfile(savepath,'UNFINISHED.lock'),'file');
  delete(fullfile(savepath,'UNFINISHED.lock'));
end;

