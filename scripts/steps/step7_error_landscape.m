% script generating error landscapes by computing the true error of reduced
% simulations vs. detailed simulations for given test parameters over differing
% basis sizes.
%
% The generated cell array 'output' can later be visualized with
% plot_error_landscape().
%
% required variables that need to be set:
%    - 'model'
%    - 'detailed_data'
%    - 'Nsize'       : number of reduced basis sizes to be tested.
%    - 'Msize'       : number of colletaral reduced basis sizes to be tested.
%    - 'csample'     : vector of increasing values between 0 and 1 indicating a
%                      set of reduced basis sizes to be tested for coupled
%                      plots.  See also couple_N_and_M_by_c() for more
%                      information.
%    - 'cdescr'      : description text for coupled plots.
%    - 'mu_set_size' : default size of test parameter set.
%    - 'model.RB_detailed_train_savepath' : path to detailed simulations in
%                       for train parameter set.
%    - 'params.step7_outputfile' : file name where the generated output
%                       structure is saved.
%
%  optional variables that can be set:
%    - 'error_plots' : cell array of strings specifying which test cases shall
%                      be handled by the script. Possible choices are
%       -# 'train_set' -- error landscape over N and M for training parameter
%                       set
%       -# 'train_set_coupled' -- error plot over c where c is the parameter
%                       of couple_N_and_M_by_c() for training parameter set
%       -# 'error' -- error landscape over N and M
%       -# 'error_coupled -- error plot over c where c is the parameter of
%                       couple_N_and_M_by_c()
%       -# 'error_ei' -- error landscape over N and M for error between
%                       reduced and detailed_simulation with empirically
%                       interpolated operators.
%       -# 'error_ei_coupled' -- error landscape over c where c is the
%                       parameter of couple_N_and_M_by_c() for error between
%                       reduced and detailed_simulation with empirically
%                       interpolated operators.
%       -# 'error_ei_rb_proj' -- error landscape over N and M for error between
%                       reduced and detailed_simulation with empirically
%                       interpolated operators with a subsequent projection on
%                       the reduced basis space in each time step.
%       -# 'error_ei_rb_proj_coupled' -- error landscape over c where c is the
%                       parameter of couple_N_and_M_by_c() for error between
%                       reduced and detailed_simulation with empirically
%                       interpolated operators with a subsequent projection on
%                       the reduced basis space in each time step.
%       .
%                    By default all error landscapes are selected.


if ~exist('error_plots','var')
  error_plots = {};
end

disp('warning: takes a few hours!');

if model.verbose < 20 && ~model.debug
  warning('off','MATLAB:nearlySingularMatrix');
end

sample_size = [Nsize, Msize];

% produce reduced data
reduced_data = gen_reduced_data(model,detailed_data);
i=0;
sconf = struct('pf',[],'max',[],'sample_size',[],'mode',[],'testdir',[],...
               'run_name_infix',[],'pf_descr',[],'testdir_infix',[]);

if isempty(error_plots) || ismember('train_set', error_plots) ...
    || ismember('train_set_coupled', error_plots)
  if ~exist(fullfile(rbmatlabtemp, model.RB_detailed_train_savepath), 'dir') ...
      || exist(fullfile(rbmatlabtemp, model.RB_detailed_train_savepath, 'UNFINISHED.lock'), 'file')
    save_detailed_simulations(model, detailed_data, ...
                              detailed_data.RB_info.M_train, ...
                              model.RB_detailed_train_savepath);
  end
end

% disable error estimator
old_model.enable_error_estimator = model.enable_error_estimator;
old_model.Mstrich                = model.Mstrich;
model.Mstrich                    = 0;
model.enable_error_estimator     = 0;

if isempty(error_plots) || ismember('train_set', error_plots)
  % sconf = struct('pf',repmat({[]}, 6,1));
  % error landscape over N and M for training set
  i=i+1;
  sconf(i).pf             = { 'N', 'M' };
  sconf(i).max            = [model.Nmax, model.Mmax];
  sconf(i).sample_size    = sample_size;
  sconf(i).mode           = 'error';
  sconf(i).testdir        = model.RB_detailed_train_savepath;
  sconf(i).run_name_infix = '';
end

if isempty(error_plots) || ismember('train_set_coupled', error_plots)
  % error landscape over c for training set, where c couples N and M by
  % (N,M)=c*(Nmax,Mmax)
  i=i+1;
  sconf(i).pf             = { @couple_N_and_M_by_c };
  sconf(i).samples        = {csample};
  sconf(i).pf_descr       = {cdescr};
  sconf(i).mode           = 'error';
  sconf(i).testdir        = model.RB_detailed_train_savepath;
  sconf(i).run_name_infix = 'coupled_';
end

if isempty(error_plots) || ismember('error', error_plots)
  % error landscape over N and M for a new test set of mu vectors.
  i=i+1;
  sconf(i).pf             = { 'N', 'M' };
  sconf(i).max            = [model.Nmax, model.Mmax];
  sconf(i).sample_size    = sample_size;
  sconf(i).mu_set_size    = mu_set_size;
  sconf(i).mode           = 'error';
  sconf(i).testdir_infix  = num2str(sconf(i).mu_set_size);
  sconf(i).run_name_infix = '';
end

if isempty(error_plots) || ismember('error_coupled', error_plots)
  % error landscape over c for a test set of mu vectors, where c couples N
  % and M by (N,M)=c*(Nmax,Mmax)
  i=i+1;
  sconf(i).pf             = { @couple_N_and_M_by_c };
  sconf(i).samples        = {csample};
  sconf(i).mu_set_size    = mu_set_size;
  sconf(i).pf_descr       = {cdescr};
  sconf(i).mode           = 'error';
  sconf(i).testdir_infix  = num2str(sconf(i).mu_set_size);
  sconf(i).run_name_infix = 'coupled_';
end

if isempty(error_plots) || ismember('error_ei', error_plots)
  % error landscape for error between reduced and detailed simulation with
  % empirically interpolated operators over N and M for a new test set of
  % mu vectors.
  i=i+1;
  sconf(i).pf             = { 'N', 'M' };
  sconf(i).max            = [model.Nmax, model.Mmax];
  sconf(i).sample_size    = sample_size;
  sconf(i).mu_set_size    = mu_set_size;
  sconf(i).mode           = 'error_to_ei';
  sconf(i).testdir_infix  = ['ei_', num2str(sconf(i).mu_set_size)];
  sconf(i).run_name_infix = '';
end

if isempty(error_plots) || ismember('error_ei_coupled', error_plots)
  % error landscape for error between reduced and detailed simulation with
  % empirically interpolated operators over c for a test set of mu vectors,
  % where c couples N and M by (N,M)=c*(Nmax,Mmax)
  i=i+1;
  sconf(i).pf             = { @couple_N_and_M_by_c };
  sconf(i).samples        = {csample};
  sconf(i).pf_descr       = {cdescr};
  sconf(i).mu_set_size    = mu_set_size;
  sconf(i).mode           = 'error_to_ei';
  sconf(i).testdir_infix  = ['ei_', num2str(sconf(i).mu_set_size)];
  sconf(i).run_name_infix = 'coupled_';
end

if isempty(error_plots) || ismember('error_ei_rb_proj', error_plots)
  % error landscape for error between reduced and detailed simulation with
  % empirically interpolated operators and an additional projection on the
  % RB space after each time step, over N and M for a new test set of mu
  % vectors.
  i=i+1;
  sconf(i).pf             = { 'N', 'M' };
  sconf(i).max            = [model.Nmax, model.Mmax];
  sconf(i).sample_size    = sample_size;
  sconf(i).mu_set_size    = mu_set_size;
  sconf(i).mode           = 'error_to_ei_rb_proj';
  sconf(i).testdir_infix  = ['ei_rb_proj_', num2str(sconf(i).mu_set_size)];
  sconf(i).run_name_infix = '';
end

if isempty(error_plots) || ismember('error_ei_rb_proj_error_plots', error_plots)
  % error landscape for error between reduced and detailed simulation with
  % empirically interpolated operators and an additional projection on the
  % RB space after each time step, over c for a test set of mu vectors,
  % where c couples N and M by (N,M)=c*(Nmax,Mmax)
  i=i+1;
  sconf(i).pf             = { @couple_N_and_M_by_c };
  sconf(i).samples        = {csample};
  sconf(i).mu_set_size    = mu_set_size;
  sconf(i).pf_descr       = {cdescr};
  sconf(i).mode           = 'error_to_ei_rb_proj';
  sconf(i).testdir_infix  = ['ei_rb_proj_', num2str(sconf(i).mu_set_size)];
  sconf(i).run_name_infix = 'coupled_';
end

output = cell(size(sconf));
% produce error landscapes for all sconf configurations above
parfor i = 1:length(sconf)
  tmp_params = params;
  range_params = [];
  range_params.plot_fields = sconf(i).pf;
  if isempty(sconf(i).max)
    range_params.samples = sconf(i).samples;
  else
    range_params.max         = sconf(i).max;
    range_params.sample_size = sconf(i).sample_size;
  end
  if ~isempty(sconf(i).pf_descr)
    range_params.plot_field_descr = sconf(i).pf_descr;
  end
  if ~isempty(sconf(i).mu_set_size)
    range_params.mu_set_size = sconf(i).mu_set_size;
  end
  if ~isempty(sconf(i).testdir)
    testdir = sconf(i).testdir;
  else
    testdir = [ tmp_params.model_type,'_',sconf(i).testdir_infix,'_test' ];
  end
  tmp_params.run_name    = [testdir, '_', ...
                        sconf(i).run_name_infix, 'step7'];
  tmp_params.mode = sconf(i).mode;
  % tmp_params.mode = 'check';

  output{i}=stochastic_error_estimation(model, detailed_data, ...
                                        reduced_data, testdir, ...
                                        range_params, tmp_params);
%  save([tmp_params.step7_outputfile,'_tmp'], 'output');
end

model.enable_error_estimator = old_model.enable_error_estimator;
model.Mstrich                = old_model.Mstrich;
clear('old_model');

save(params.step7_outputfile, 'output');

