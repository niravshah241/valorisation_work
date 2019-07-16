function [detailed_data, reduced_data] = RB_greedy_extension(model,detailed_data)
%function [detailed_data,reduced_data] = ...
%             RB_greedy_extension(model,detailed_data);
%  function performing a greedy search loop for reduced basis
%  generation.
%
%  In each iteration the parameter vector with maximum posterior estimator or
%  true error is selected for "refinement". An orthonormal basis is generated.
%
%  The detailed_data is expected to be reasonably filled, that rb_simulations
%  are possible with the given arguments, i.e.  gen_reduced_data(model,
%  detailed_data). That means, that it also can contain an existing basis,
%  which is then extended in this loop.  the offline_data for after the last
%  extension is optionally returned.
%
%  The stopping of the extention is performed, if
%
%    - 'model.RB_stop_epsilon' is reached: basis vectors are added until
%             the maximum posterior estimator is smaller than this
%             value. M_val can be ommited in this case.
%    - 'model.RB_stop_max_val_train_ratio': basis vectors are added until the
%             ratio r of the max estimator on a validation set M_val with the
%             estimator on the training set M is larger than r_max. M_val is
%             required in this case
%    - 'model.RB_stop_timeout': basis vectors are added until toc gives a
%             time exceeding a maximum level given
%    - 'model.RB_stop_Nmax': If given number of basis vectors in full
%             basis is obtained.
%    - 'model.p_part_early_refinement' is set to one and the estimated
%             number of basis vectors necessary to reach etol is higher
%             than 'RB_stop_Nmax'.
%
% The following fields are created/extended in 'detailed_data.RB_info':
%    - 'M_first_errs' : a vector with the initial error-estimator for each
%               M_train. Is determined to check, whether any
%               extension has to be performed or not.
%    - 'stopped_on_epsilon' : flag indicating whether epsilon-stop-crit is reached
%    - 'stopped_on_timeout': flag indicating whether timeout is reached
%    - 'stopped_on_Nmax': flag indicating whether Nmax is reached
%    - 'stopped_on_max_val_train_ratio': flag indicating whether
%               epsilon-stop-crit is reached
%    - 'stopped_on_empty_extension' : flag indicating, whether the last
%               extension failed due to accuracy problems, i.e. new
%               basis-vector already in span of existing.
%    - 'stopped_on_Nlimit_estimation': flag indicating ifthe extension
%               process was stopped, because it was estimated that more
%               than Nmax basis vectoras are needed to reach etol.
%               Partition of the parameter domain can be initiated.
%    - 'max_err_sequence' : the decreasing maximum posterior estimators or true
%               errors in each iteration, which are present AFTER the
%               respective basis extension. So in case of 'stopped_on_epsilon',
%               the last entry of this field should be smaller than the desired
%               epsilon.
%    - 'mu_sequence' : sequence of mu_values chosen for
%               basis extension i.e. list of worst-estimator mu_values as the
%               columns of the matrix mu_values.
%    - 'r_value_sequence' : ratios of val/train error-estimator in case of
%               validation set providing, AFTER the respective basis-extension.
%    - 'toc_value_sequence' : toc_values after each basis extension
%    - 'M_last_errs' : the posterior estimators for all mu in M_train at final
%               time AFTER the last respective basis extension
%
%  Required fields in 'detailed_data.RB_info':
%     - 'M_train' : matrix with parameter-vectors as columns, over which greedy
%                   search is to be performed.
%     - 'M_val'   : matrix with parameter-vectors as columns, which are used
%                   for validation. Only required, if the field
%                   'model.stop_on_max_val_train_ratio' is set)
%
%  The following fields of detailed_data are created/extended:
%     - RB      : is extended by the number of basis-vectors, overall being the
%                 orthonormal basis
%     - RB_info : the fields of this info-structure are extended, as described
%                 above
%
% required fields of model:
%    RB_extension_algorithm: function pointer to either
%                         RB_extension_max_error_snapshot() or
%                         RB_extension_PCA_fixspace()
%    RB_error_indicator: either 'estimator' or 'error' indicating,
%                  whether an aposteriori indicator or true error
%                  is to be used
%    RB_detailed_val_savepath : path to the detailed data of the
%                  validation set in validation mode. Is generated
%                  if not available.
%    RB_detailed_train_savepath : path to the detailed data of the
%                  training set. Is generated if not available. Is
%                  only required in case of 'error' mode instead of
%                  'estimator'
%
% additional fields as required by the selected
% extension algorithm and the reduced-basis-simulation must be set in model.

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


% Bernard Haasdonk 26.5.2007

if ~isfield(model,'RB_stop_timeout')
  model.RB_stop_timeout             = inf;
end;

if ~isfield(model,'RB_stop_Nmax')
  model.RB_stop_Nmax                = inf;
end;

if ~isfield(model,'RB_stop_epsilon')
  model.RB_stop_epsilon             = 0;
end;

if ~isfield(model,'RB_stop_max_val_train_ratio')
  model.RB_stop_max_val_train_ratio = inf;
end;

% check, whether val-error is wanted
check_val_error = 0;
if model.RB_stop_max_val_train_ratio < inf
  check_val_error                   = 1;
end;

if isfield(model,'p_part_early_refinement')
  enable_p_partition                = 1;
else
  enable_p_partition                = 0;
end

% make sure, that detailed-simulations for validation set are available
% in case of active validation criterion and wanted true error 
if check_val_error
  if isempty(detailed_data.RB_info.M_val)
    error ('M_val must be provided in case of max_val_train_ratio mode!');
  end;
  if isequal(model.RB_error_indicator,'error')
    detailed_val_savepath = model.RB_detailed_val_savepath;
    save_detailed_simulations_ptr = model.save_detailed_simulations;
    save_detailed_simulations_ptr(model,detailed_data, ...
                                  detailed_data.RB_info.M_val, ...
                                  detailed_val_savepath);
  else
    detailed_val_savepath = [];
  end;
end;

% make sure, that detailed-simulations for the training set are available in
% case of wanted true error 
if isequal(model.RB_error_indicator,'error')
  save_detailed_simulations_ptr = model.save_detailed_simulations;
  save_detailed_simulations_ptr(model,detailed_data,...
                                detailed_data.RB_info.M_train, ...
                                model.RB_detailed_train_savepath);
  detailed_train_savepath = model.RB_detailed_train_savepath;
else
  detailed_train_savepath = 'dummy-string';
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% determine mu-worst error-indicator for current training set %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

reduced_data = gen_reduced_data(model,detailed_data);

model.N = reduced_data.N;
post_errs = rb_test_indicator(model, ...
                              detailed_data, reduced_data,...
                              detailed_data.RB_info.M_train,...
                              detailed_train_savepath);
% store real error as well, if estimator and error are computed for
% 'ei_estimator_test'
if isequal(model.RB_error_indicator,'ei_estimator_test')
  post_errs_sequence = post_errs.error;
  post_est_sequence  = post_errs.estimator;
  post_errs          = post_errs.estimator;
end

[max_errs, mu_inds] = max(post_errs);
max_err             = max_errs(1);      % maximum error (estimation)
mu_ind              = mu_inds(1);       % parameter index in RB_info.M_train

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize loop-break conditions %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
stopped_on_epsilon             = 0;
stopped_on_timeout             = 0;
stopped_on_Nmax                = 0;
stopped_on_max_val_train_ratio = 0;
stopped_on_empty_extension     = 0;
stopped_on_Nlimit_estimation   = 0;
N_limit_est_av=[]; % average over all estimations (pessimistic)
continue_loop = 1;
if isnan(max_err)
    keyboard
end;

% check whether epsilon has been reached
if max_err < model.RB_stop_epsilon
  continue_loop      = 0;
  stopped_on_epsilon = 1;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize output quantities %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
max_pe     = [];
mu_values  = [];
mu_ind_seq = [];
toc_values = [];
if check_val_error
  r_values = [];
end;
N_limit_est = [];% estimated amount of N needed


% initialize detailed_data.RB_info
if ~isfield(detailed_data,'RB_info')
  detailed_data.RB_info                    = [];
end;

N = model.get_rb_size(model,detailed_data);

if ~isfield(detailed_data.RB_info,'max_err_sequence')
  detailed_data.RB_info.max_err_sequence   = nan(1,N);
end
if ~isfield(detailed_data.RB_info,'mu_sequence')
  detailed_data.RB_info.mu_sequence        = ...
      nan(length(model.mu_names), N);
end
if ~isfield(detailed_data.RB_info,'mu_ind_seq')
  detailed_data.RB_info.mu_ind_seq         = ...
      nan(1, N);
end

if check_val_error && ~isfield(detailed_data.RB_info,'r_value_sequence')
  detailed_data.RB_info.r_value_sequence   = nan(1,N);
end

if ~isfield(detailed_data.RB_info,'toc_value_sequence')
  detailed_data.RB_info.toc_value_sequence = nan(1,N);
end

detailed_data.RB_info.M_first_errs = post_errs;

% perform extension loop:
while continue_loop
  % Assume that mu_ind points to the current update mu-vector
  disp(['Detected maximum error prediction ',num2str(max_err),...
        ' for mu=[',...
        num2str(detailed_data.RB_info.M_train(:,mu_ind)'),']']);

  % Use original params with new mu such that filecaching works in the
  % extension routine!
  mu = detailed_data.RB_info.M_train(:,mu_ind);

  oldmu        = model.get_mu(model);
  model        = model.set_mu(model, mu);
  [RBext, par] = model.RB_extension_algorithm(model, ...
                                              detailed_data);
  model        = model.set_mu(model, oldmu);

  % if empty: no further processing has to be performed
  if isempty(RBext)
    stopped_on_empty_extension = 1;
    continue_loop              = 0;
  else
    stopped_on_empty_extension = 0;
    if model.verbose > 29
      disp(par);
    end;
    max_pe     = [max_pe, max_err];
    mu_ind_seq = [mu_ind_seq, mu_ind];
    mu_values  = [mu_values, detailed_data.RB_info.M_train(:,mu_ind)];
    toc_values = [toc_values, toc];
    RB = [model.get_rb_from_detailed_data(detailed_data), RBext];
    % perform permutation that maintains original order if possible:
    %RBo = model_orthonormalize_qr(model,detailed_data,RB,eps);
    %K = model.inner_product(model,detailed_data,RBo,RB);
    %[m,i] = max(abs(K));
    %RB = RBo(:,i);
%    if model.debug
    %  K  = model.inner_product(model,detailed_data,RB,RB);
    %  if max(max(abs(K-eye(size(K)))))>1e-10
    %    	%error('error in orthonormal basis extension!');
    %  end;
%   % end;
    detailed_data = model.set_rb_in_detailed_data(...
        detailed_data,...
        RB);
    %detailed_data.RB = [detailed_data.RB, RBext];
    N = model.get_rb_size(model,detailed_data);
    detailed_data.N = N;
    disp(['Extended RB to length ', num2str(N)]);

  end;


  %orthogonality check

  disp(['checking othogonality with tolerance',num2str(1e-12),'...',num2str(max(max(detailed_data.RB'*detailed_data.W*detailed_data.RB-eye(size(detailed_data.RB,2)))))]);
  if(max(max(detailed_data.RB'*detailed_data.W*detailed_data.RB-eye(size(detailed_data.RB,2))))>1e-12)
    disp('orthogonalising again...')
    detailed_data.RB = model_orthonormalize_qr(model, detailed_data,detailed_data.RB);
    disp(['new orthogonality value:  ',num2str(max(max(detailed_data.RB'*detailed_data.W*detailed_data.RB-eye(size(detailed_data.RB,2)))))]);
  end


%  keyboard;

  % check on accuracy problem (non-decreasing error)
%  if length(max_pe)>1 & max_pe(end)>max_pe(end-1)
%    stopped_on_error_increase = 1;
%    continue_loop = 0;
%    disp('stopped on error-increase! Please check!');
%    keyboard;
%  else
%    stopped_on_error_increase = 0;
%  end;


  if ~stopped_on_empty_extension
    % determine error estimators for extended basis for next loop
    reduced_data = gen_reduced_data(model,detailed_data);
    model.Nmax   = model.get_rb_size(model,detailed_data);
    model.N      = model.Nmax;
    %    model = set_mu(detailed_data.RB_info.M_train(:,mu_ind), ...
    %                   model);
    post_errs = rb_test_indicator(model,...
                                  detailed_data, reduced_data,...
                                  detailed_data.RB_info.M_train,...
                                  detailed_train_savepath);
    if isequal(model.RB_error_indicator,'ei_estimator_test')
      post_errs_sequence = [post_errs_sequence, post_errs.error];
      post_est_sequence  = [post_est_sequence, post_errs.estimator];
      post_errs          = post_errs.estimator;
    end
    [max_errs, mu_inds] = max(post_errs);
    max_err             = max_errs(1);
    mu_ind              = mu_inds(1);

%    model.Nmax = size(detailed_data.RB,2);
%    model.N = model.Nmax;
%    [mu_ind, post_errs] = detect_mu_worst(detailed_data,model);
%    max_err = max(post_errs);

    % check loop-break conditions
    if max_err < model.RB_stop_epsilon
      continue_loop      = 0;
      stopped_on_epsilon = 1;
    end;

    if model.get_rb_size(model,detailed_data) >= model.RB_stop_Nmax
      continue_loop   = 0;
      stopped_on_Nmax = 1;
    end;

    if check_val_error
      val_est = rb_test_indicator(model,...
                                  detailed_data,reduced_data,...
                                  detailed_data.RB_info.M_val, ...
                                  detailed_val_savepath);

      r_values = [r_values, max(val_est) / max_err];
      %debug
      disp(['val_train ratio:  ',num2str(r_values(end))]);
      %end debug
      if (r_values(end) > model.RB_stop_max_val_train_ratio)
        stopped_on_max_val_train_ratio = 1;
        continue_loop = 0;
      end;
    end;

    %in models using adaptive p-partition early stopping:
    %an estimater extrapolates the curve N->epsilon and estimates if the
    %epsilon_tol can be reached with a basis size smaller than Nmax. If no,
    %the field RB_stopped_on_Nlimit_estimation is set to 1


    if enable_p_partition && model.p_part_early_refinement
      %estimate N total needed to reach etol:
      err_seq=max_pe;
      etol = model.RB_stop_epsilon;
      Nmax = model.RB_stop_Nmax;

      s=20; %stepsize backwards


      if(length(err_seq)>(s))

        %for i=(s+1):length(err_seq);
        i=length(err_seq);
        a=(log(err_seq(i-s)/err_seq(i)))/(s);
        N_limit_est=[N_limit_est,i-1/a*log(etol/err_seq(i))];
        %end;



        %for i=1:length(N_limit_est);
        %summe = summe+N_limit_est(i);
        if(length(N_limit_est_av>0))
          N_1=length(N_limit_est_av);
          newAv=N_1/(N_1+1)*N_limit_est_av(end)+(1/(N_1+1))*N_limit_est(end);
          N_limit_est_av=[N_limit_est_av,newAv];

        else
          N_limit_est_av=N_limit_est(end);
        end
        %end;
        disp(['actual estimation of N needed to reach etol: ',num2str(N_limit_est_av(end))]);
        if(N_limit_est_av(end)>Nmax)
          stopped_on_Nlimit_estimation = 1;
          continue_loop                = 0;
        end
      end


    end %of P-Part_early_refinement

  end; % of ~stopped_on_empty_extension

  if toc > model.RB_stop_timeout
    stopped_on_timeout = 1;
    continue_loop      = 0;
  end;

end; % end while loop

% prepare output data : RB is already updated 

RB_info                                = detailed_data.RB_info;

RB_info.stopped_on_epsilon             = stopped_on_epsilon;
RB_info.stopped_on_max_val_train_ratio = stopped_on_max_val_train_ratio;
RB_info.stopped_on_timeout             = stopped_on_timeout;
RB_info.stopped_on_Nmax                = stopped_on_Nmax;
RB_info.stopped_on_empty_extension     = stopped_on_empty_extension;
RB_info.stopped_on_Nlimit_estimation   = stopped_on_Nlimit_estimation;

RB_info.M_last_errs      = post_errs;
RB_info.max_err_sequence = [RB_info.max_err_sequence, max_pe];
RB_info.mu_ind_seq       = [RB_info.mu_ind_seq, mu_ind_seq];
RB_info.mu_sequence      = [RB_info.mu_sequence, mu_values];
if isequal(model.RB_error_indicator,'ei_estimator_test')
  RB_info.post_errs_sequence = post_errs_sequence;
  RB_info.post_est_sequence  = post_est_sequence;
end
if check_val_error
  RB_info.r_value_sequence = [RB_info.r_value_sequence, r_values];
end;
RB_info.toc_value_sequence = [RB_info.toc_value_sequence, toc_values];

if enable_p_partition && model.p_part_early_refinement
   RB_info.N_limit_est    = N_limit_est;
   RB_info.N_limit_est_av = N_limit_est_av;
end


detailed_data.RB_info = RB_info;

