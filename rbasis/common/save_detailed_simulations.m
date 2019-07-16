function save_detailed_simulations(model,model_data,M,savepath)
% function save_detailed_simulations(model,model_data,M,savepaths)
% perform loop over detailed simulations and save results or check
% consistency with existing saved results.
%
% -# If the given path does not exist, it is generated and a set of
% detailed simulations for all parameters mu (columns of M).
% is computed and the results are stored in the path specified in
% 'savepath'.
% Files to be generated in the target path are:
%     - 'settings.mat'  :containing 'model', 'grid', 'M', 'path', ... such that
%                        the directory can be regenerated with this file and
%                        the current function
%     - 'detail1.mat'   :containing the simulation data for first `mu`
%                        vector in 'M'
%     - 'detail1.mat'   :containing the simulation data for second `mu`
%                        vector in 'M'
%     - ....
%     - 'detail123.mat' : ...
%     .
% -# If the path exists, a check on consistency is performed by checking the
% 'model' and 'Mtrain' fields and checking whether the previous computation has
% been finished. If it is not, and the current settings are consistent with the
% old ones, the computations are resumed.
%
% optional fields of model:
%   force_delete : if this option is set to true (the default), directories
%                  with data inconsistent with the current settings are
%                  deleted.
%
% Parameters:
%  M       : is a matrix of `mu` vectors (column vectors) for which detailed
%            simulations shall be computed and stored.
%  savepath: is the directory name relative to 'RBMATLABTEMP' where the
%            detailed simulations shall be stored.
%

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


% Bernard Haasdonk 29.3.2007

sp = fullfile(rbmatlabtemp,savepath);

if ~isfield(model, 'force_delete')
  model.force_delete = 1;
end

% generate directory if not available
if ~exist(sp,'dir')
  disp('data directory does not exist, creating .... ');
  [p,s,e] = fileparts(sp);
  mkdir(p,[s,e]);
  % now it must exist
  if ~exist(sp,'dir')
    error('error in generating savepath!');
  end;
elseif exist(fullfile(sp, 'settings.mat'), 'file')
  % check consistency of file
  tmp = load(fullfile(sp,'settings.mat'));
  % fields to skip in comparison
  ignorelist = {'Mmax','MM','Mstrich','ei_Mmax','k','Mmax','M','N','Nmax',...
              'RB_error_indicator','RB_stop_Nmax',...
              'RB_stop_timeout','RB_detailed_test_savepath',...
              'RB_generation_mode','RB_train_rand_seed','RB_train_size',...
              'test_N_samples','RB_stop_epsilon','RB_M_val_size',...
              'RB_numintervals','RB_refinement_mode',...
              'RB_refinement_theta','RB_max_refinement_level',...
              'RB_stop_max_val_train_ratio','RB_val_rand_seed',...
                'ei_stop_on_Mmax','CRB_generation_mode',...
                'extend_crb', ... 
                'CRB_basis_filename','ei_space_operators',...
              'Msamples','ei_operator_savepath','ei_target_error',...
                'ei_numintervals', 'ei_detailed_savepath','verbose','debug', ...
               'velocity_coefficients_ptr','enable_error_estimator',...
               'force_delete', 'ei_time_splits', 'adaptive_time_split_mode', ...
               'ei_stop_epsilon', 'ei_stop_epsilon_first', 'minimum_time_slice',  ...
       'time_split_Mmax'};
  ignorelist = [ignorelist, model.mu_names,...
                model.filecache_ignore_fields_in_model];

  [iseq, a,b,c] = structcmp(model,tmp.model,ignorelist);
  if ~iseq
    disp('fields of model and stored model differs:')
    disp('differing fields:');
    cellfun(@disp,a);
    disp('additional fields in model:');
    cellfun(@disp,b);
    disp('additional fields in stored model:');
    cellfun(@disp,c);
    if ~model.debug && model.force_delete
      disp(['parameters in precomputed data and current are inconsistent!',...
            ' I am deleting the old data now!']);
      delete(fullfile(sp,'*.mat'))
      ls(sp);
    else
      error(['parameters of precomputed data and current',...
             ' simulation are inconsistent! Please delete files in path ', ...
              sp, ' and restart!']);
    end
  elseif ~isequal(M,tmp.M)
    if ~model.debug && model.force_delete
      disp(['M of precomputed and current data is inconsistent!',...
            ' I am deleting the old data now!']);
      delete(fullfile(sp,'*.mat'))
      ls(sp);
    else
      error(['M of precomputed data and current',...
             ' simulation are inconsistent! Please delete and restart!']);
    end
  else
    if exist(fullfile(sp,'UNFINISHED.lock'),'file')
      if ~model.debug
        disp(['detected UNFINISHED.lock in target directory, ',...
              'resuming previous computations...']);
      else
        error(['Previous computation of detailed simulations seems not', ...
          ' to be finished!!! please delete and restart!', ...
          ' savepath=', sp]);
      end
    else
      disp('skipping detailed computations and using stored results');
      return;
    end
  end
else
  disp(['stray directory found! I am cleaning this one, now!']);
  delete(fullfile(sp,'*.mat'))
  ls(sp);
end;

save(fullfile(sp,'UNFINISHED.lock'),'model');

% save general settings
save(fullfile(sp,'settings.mat'),'M','model','savepath');

if ~isfield(model, 'num_cpus')
  model.num_cpus = 4;
end
num_cpus = model.num_cpus;

% save simulation_data
npar = size(M,2);
for mout = 0:floor((npar-1)/num_cpus)
  sim_data_c = cell(1, npar);
  tictoc_c   = cell(1, npar);
  parfor mu = (mout*num_cpus+1):min(npar,(mout+1)*num_cpus)

    filename = fullfile(sp,['detail',num2str(mu),'.mat']);
    if ~exist(filename, 'file')
      tmp_model = model;
      disp(['processing parameter vector ',num2str(mu),'/',num2str(npar)]);
      tmp_model = tmp_model.set_mu( model,M(:,mu) );
      tic;
      sim_data_c{mu} = detailed_simulation(tmp_model, model_data);
      tictoc_c{mu} = toc;
    else
      disp(['file exists: skipping parameter vector ', ...
        num2str(mu), '/', num2str(npar)]);
    end
  end;
  for mu=(mout*num_cpus+1):min(npar,(mout+1)*num_cpus)

    sim_data = sim_data_c{mu};
    tictoc   = tictoc_c{mu};
    filename = fullfile(sp,['detail',num2str(mu),'.mat']);
    if ~exist(filename, 'file')
      save(filename,'sim_data','tictoc');
    end
  end
end

delete(fullfile(sp,'UNFINISHED.lock'));


