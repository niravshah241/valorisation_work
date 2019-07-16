% script generating landscapes plot data by computing the error estimator of
% reduced simulations for given test parameters over differing basis sizes.
%
% The generated cell array 'output' can later be visualized with
% plot_error_landscape().
%
% required variables that need to be set:
%    - 'model'
%    - 'detailed_data'
%    - 'Nsize'       : number of reduced basis sizes to be tested.
%    - 'Msize'       : number of colletaral reduced basis sizes to be tested.
%    - 'Mstrich_samples' : vector of values that are assigned to M\' for plots
%                      where M\' varies.
%    - 'csample'     : vector of increasing values between 0 and 1 indicating a
%                      set of reduced basis sizes to be tested for coupled
%                      plots.  See also couple_N_and_M_by_c() for more
%                      information.
%    - 'cdescr'      : description text for coupled plots.
%    - 'mu_set_size' : default size of test parameter set.
%    - 'model.RB_detailed_train_savepath' : path to detailed simulations
%                       for train parameter set.
%    - 'params.step8_outputfile' : file name where the generated output
%                       structure is saved.
%
%  optional variables that can be set:
%    - 'estimator_plots' : cell array of strings specifying which test cases
%                          shall be handled by the script. Possible choices are
%       -# 'Mstrich_1_coupled' -- plot of error estimates over 'c', where 'c'
%                          couples 'N' and 'M' by couple_N_and_M_by_c() and
%                          'Mstrich' is set to 1 .
%       -# 'Mstrich_max_coupled' -- plot of error estimates over 'c', where 'c'
%                          couples 'N' and 'M' by couple_N_and_M_by_c() and
%                          'Mstrich' is set to maximum value in
%                          'Mstrich_samples'.
%       -# 'Mstrich_1' -- landscape plot over 'N' and 'M' and 'Mstrich' is
%                          set to 1.
%       -# 'Mstrich_max' -- landscape plot over 'N' and 'M' and 'Mstrich' is
%                          set to maximum value in array 'Mstrich_samples'.
%       -# 'Mmax_train' -- landscape plot over 'N' and 'Mstrich' with 'M' set to
%                          maximum possible value for test parameter set used
%                          in reduced basis generation.
%       -# 'coupled' -- landscape plot over 'c' and 'Mstrich' where 'c' couples
%                          'N' and 'M' by couple_N_and_M_by_c().
%       -# 'Mmax' -- landscape plot over 'N' and 'Mstrich' with 'M' set to
%                          maximum possible value.
%       -# 'Nmax' -- landscape plot over 'M' and 'Mstrich' with 'N' set to
%                          maximum reduced basis size.
%       .
%
%                    By default all estimator landscapes are selected.

if ~exist('estimator_plots','var')
  estimator_plots = {};
end

if ~isfield(params, 'N')
  params.N = model.Nmax;
end

disp('warning: takes a few hours!');

if model.verbose < 20 && ~model.debug
  warning('off','MATLAB:nearlySingularMatrix');
end

% reduced basis sizes to be tested (depends on Nsize)
Nsamples        = 1+round((0:Nsize-1) .* (model.Nmax-1) ./ (Nsize-1));
% collateral reduced basis sizes to be tested (depends on Msize)
Msamples        = 1+round((0:Msize-1) .* (model.Mmax-1) ./ (Msize-1));

% produce reduced data
reduced_data = gen_reduced_data(model,detailed_data);

i=0;
sconf = struct('pf',[],'samples',[],'pf_descr',[],'Mstrich',[],...
               'M',[],'N',[],'testdir_infix',[],'mu_set',[],...
               'max',[]);
if isempty(estimator_plots) || ismember('Mstrich_1_coupled', estimator_plots)
  % plot of error estimates over 'c', where 'c' couples 'N' and 'M' by `(N,M) =
  % c (Nmax, Mmax)` and 'Mstrich' is set to 1
  i=i+1;
  sconf(i).pf            = { @couple_N_and_M_by_c };
  sconf(i).samples       = {csample};
  sconf(i).pf_descr      = {cdescr};
  sconf(i).Mstrich       = 1;
  sconf(i).testdir_infix = 'coupled_Mstrich_1_';
end

if isempty(estimator_plots) || ismember('Mstrich_max_coupled', estimator_plots)
  % plot of error estimates over 'c', where 'c' couples 'N' and 'M' by `(N,M) =
  % c (Nmax, Mmax)` and 'Mstrich' is set to its maximum value
  i=i+1;
  sconf(i).pf            = { @couple_N_and_M_by_c };
  sconf(i).samples       = {csample};
  sconf(i).pf_descr      = {cdescr};
  sconf(i).Mstrich       = max(Mstrich_samples);
  sconf(i).testdir_infix = 'coupled_Mstrich_Mmax_';
end

if isempty(estimator_plots) || ismember('Mmax_train', estimator_plots)

  if ~exist(fullfile(rbmatlabtemp, model.RB_detailed_train_savepath), 'dir') ...
    || exist(fullfile(rbmatlabtemp, model.RB_detailed_train_savepath, 'UNFINISHED.lock'), 'file')
    save_detailed_simulations(model, detailed_data, ...
                              detailed_data.RB_info.M_train, ...
                              model.RB_detailed_train_savepath);
  end

  tmp = load(fullfile(rbmatlabtemp, model.RB_detailed_train_savepath,'settings'));
  i=i+1;
  sconf(i).pf            = { 'N' };
  sconf(i).samples       = { Nsamples };
  sconf(i).mu_set        = tmp.M;
  sconf(i).pf_descr      = {cdescr};
  sconf(i).Mstrich       = max(Mstrich_samples);
  sconf(i).M             = model.Mmax;
  sconf(i).testdir_infix = 'N_Mstrich_Mmax_train_set';
end

if isempty(estimator_plots) || ismember('Mstrich_1', estimator_plots)
  % plot of error estimates over 'N' and 'M', where and 'Mstrich' is set to
  % one.
  i=i+1;
  sconf(i).pf            = { 'N', 'M' };
  sconf(i).max           = [ model.Nmax, model.Mmax ];
  sconf(i).sample_size   = [ Nsize, Msize ];
  sconf(i).Mstrich       = 1;
  sconf(i).testdir_infix = 'Mstrich_1_';
end

if isempty(estimator_plots) || ismember('Mstrich_max', estimator_plots)
  % plot of error estimates over 'N' and 'M', where and 'Mstrich' is set to
  % maximum.
  i=i+1;
  sconf(i).pf            = { 'N', 'M' };
  sconf(i).max           = [ model.Nmax, model.Mmax ];
  sconf(i).sample_size   = [ Nsize, Msize ];
  sconf(i).Mstrich       = max(Mstrich_samples);
  sconf(i).testdir_infix = 'Mstrich_max_';
end

if isempty(estimator_plots) || ismember('coupled', estimator_plots)
  % plot of error estimates over 'c' and 'Mstrich', where 'c' couples 'N' and
  % 'M' by `(N,M) = c (Nmax, Mmax)`.
  i=i+1;
  sconf(i).pf            = { @couple_N_and_M_by_c, 'Mstrich' };
  sconf(i).samples       = {csample, Mstrich_samples};
  sconf(i).pf_descr      = {cdescr, '\hat{M}'};
  sconf(i).testdir_infix = 'coupled_Mstrich_';
end

if isempty(estimator_plots) || ismember('Mmax', estimator_plots)
  % plot of error estimates over 'N' and 'Mstrich'
  i=i+1;
  sconf(i).pf            = { 'N', 'Mstrich' };
  sconf(i).samples       = { Nsamples, Mstrich_samples };
  sconf(i).pf_descr      = { 'N', '\hat{M}' };
  sconf(i).testdir_infix = 'N_Mstrich_';
  sconf(i).M             = model.Mmax;
end

if isempty(estimator_plots) || ismember('Nmax', estimator_plots)
  % plot of error estimates over 'M' and 'Mstrich'
  i=i+1;
  sconf(i).pf            = { 'M', 'Mstrich' };
  sconf(i).samples       = { Msamples, Mstrich_samples };
  sconf(i).pf_descr      = { 'M', '\hat{M}' };
  sconf(i).testdir_infix = 'M_Mstrich_';
  sconf(i).N             = params.N;
end

output = cell(1,length(sconf));
% produce plots specified by the configuration structure 'sconf'
for i = 1:length(sconf)
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
  if ~isempty(sconf(i).Mstrich)
    model.Mstrich = sconf(i).Mstrich;
  end
  if ~isempty(sconf(i).M)
    model.M = sconf(i).M;
  end
  if ~isempty(sconf(i).N)
    model.N = sconf(i).N;
  end
  if ~isempty(sconf(i).mu_set)
    range_params.mu_set = sconf(i).mu_set;
    testdir = [ params.model_type, sconf(i).testdir_infix, ...
                '_train' ];
  else
    range_params.mu_set_size = mu_set_size;
    testdir = [ params.model_type, sconf(i).testdir_infix, '_', ...
                num2str(range_params.mu_set_size) ];
  end
  params.run_name    = [ testdir, '_step8'];
  params.mode = 'estimator';
  %params.mode = 'check';

  output{i}=stochastic_error_estimation(model, detailed_data, ...
                                       reduced_data, testdir, ...
                                       range_params, params);

  save([params.step8_outputfile,'_tmp'], 'output');
end
save(params.step8_outputfile, 'output');
