function output=stochastic_error_estimation(model, detailed_data, ...
                                            reduced_data, testdir, ...
                                            range_params, params)
% function stochastic_error_estimation(model, detailed_data, reduced_data, ...
%                                      testdir, range_params, params)
% stochastic estimation of error between reduced and detailed simulation over a
% test sample of `\mu` vectors
%
% This function stochastically estimates the error between reduced and detailed
% simulations for given `\mu`-vectors and various reduced and collateral basis
% sizes. The results are visualized in a surface plot for problems with CRB.
%
% If required by the users, averaged time measurements for the reduced and the
% detailed simulations are computed, too.
%
% Parameters:
%  testdir:       is a directory name where the detailed simulations are stored
%  range_params:  specifies the combinations of RB and CRB sizes to be tested.
%  params:        specifies control parameters
%
% Required fields of range_params:
%  plot_fields:   is a vector of field names or function pointers with at most
%                 2 elements over which the error landscape is computed. A
%                 reasonable choice would be '{ \'N\', \'M\' }'. The error
%                 landscape is plotted over a variation of these fields
%                 specified by other 'range_params' fields. In case of
%                 function pointers, 'model' fields are updated before each
%                 new error computation via calls to these function pointers.
%                 The required synopsis for such function pointers 'fptr' is:
%                 @code
%                   function model = fptr(model, new_value)
%                 @endcode
%
% Optional fields of range_params:
%  min:         is a vector of minimum values of 'plot_fields' variables. The
%               default value is '[1,1]'
%  max:         is a vector for maximum values of 'plot_fields' variables.
%               Either this field or 'samples' need to be set.
%  sample_size: is only used, if 'max' is set. It specifies the number of
%               sample values between 'min' and 'max' for which the error is
%               computed and plotted. The default value is 'max-min'.
%  samples:     is a cell array of scalar row vectors. The scalar values are
%               'plot_field' values for which the error is computed and tested.
%  mu_set_size: number of `\mu` vectors to be tested. If this field is not set,
%               the directory given by 'testdir' needs to hold detailed
%               simulations computed by the 'save_detailed_simulations' method,
%               or the fied 'mu_set' needs to be set.
%  mu_set:      set of `\mu` vectors to be tested. If this field is not set,
%               the directory given by 'testdir' needs to hold detailed
%               simulations computed by the 'save_detailed_simulations' method,
%               or the fied 'mu_set_size' needs to be set.
%  plot_field_descr: is a cell array of description texts for the plot fields.
%               If this field is not set it is set to 'plot_fields'.
%
% Required fields of params:
%  run_name:      character string specifying this test run. The
%                 string should be unique for every parameter
%                 combination.
%
% Optional fields of params:
%  mode:          - 'error' (default)  computes the error between reduced and
%                   detailed simulations.
%                 - 'error_to_ei' computes the error between reduced and
%                   detailed simulation with empirically approximated operators.
%                 - 'error_to_ei_rb_proj' computes the error between reduced
%                   and detailed simulation with empiricaly approximated
%                   operators and subsequent projection to the RB space.
%                 - 'estimator' computes the error estimation for reduced
%                   simulations.
%                 - 'check' computes nothing, but can be used, for a
%                    quick check wether parameters are set correctly.
%  email:         email address where a progress notification can be sent to
%  mail_interval: After 'params.mail_interval' `\mu` vectors have been
%                 processed a notification to the email address given by
%                 'params.email' is sent. (Default = '10')
%  stab_limit:    stability region - If the averaged error exceeds this value,
%                 the visualization is cropped here. This happens especially
%                 for small reduced basis sizes. (Default = '1e-2')

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

output = [];

if ~isfield(params, 'mode')
  params.mode = 'error';
end

if isequal(params.mode, 'error_ei')
  model.detailed_simulation = model.detailed_ei_simulation;
end

if isequal(params.mode, 'error_ei_rb_proj')
  model.detailed_simulation = model.detailed_ei_rb_proj_simulation;
end

if isequal(params.mode, 'estimator')
  model.RB_error_indicator = 'estimator';
end

% default value (only need for mode == 'check' when no 'mu_set_size' is given)
M = rand_uniform(1, model.mu_ranges);

if isfield(range_params, 'mu_set_size') && range_params.mu_set_size > 0
  rand('state', 123456);
  M = rand_uniform(range_params.mu_set_size, model.mu_ranges);
  if ~isequal(params.mode, 'estimator') && ...
      ~isequal(params.mode, 'check')
    save_detailed_simulations(model, detailed_data, M, testdir);
  end
elseif isfield(range_params, 'mu_set') && ~isempty(range_params.mu_set)
  M = range_params.mu_set;
end

if ~isfield(range_params, 'min')
  range_params.min = [1,1];
end

if isfield(range_params, 'max')
  if isfield(range_params, 'samples')
    warning('max and samples both set in range_params!');
  end
  ssize = range_params.sample_size;
  tsamples = cell(1,length(range_params.max));
  for i = 1:length(range_params.max)
    tsamples{i} = range_params.min(i) + ...
      round((0:ssize(i)-1) .* ...
            (range_params.max(i)-range_params.min(i))./(ssize(i)-1));
  end
else
  if ~isfield(range_params, 'samples')
    error('Either max or samples need to be set in range_params!');
  else
    if ~iscell(range_params.samples)
      error('samples need to be a cell array.');
    end
    tsamples = range_params.samples;
  end
end

pf_descr = cell(1, length(tsamples));
for i = 1:length(tsamples)
 if isfield(range_params, 'plot_field_descr') && ...
     ~isempty(range_params.plot_field_descr)
   pf_descr{i} = range_params.plot_field_descr{i};
 else
   pf_descr{i} = range_params.plot_fields{i};
 end
end

if isfield(params, 'mail_interval')
  mail_interval = params.mail_interval;
else
  mail_interval = 10;
end

% define stability-region as error being smaller than
% sqrt(diffmax * area), e.g. diffmax = 4, area = 2e-7
if isfield(params, 'stab_limit')
  stab_limit = params.stab_limit;
else
  stab_limit = 1e-2;
end

savefile = [params.run_name];
tmpfile  = [savefile,'_tmp.mat'];
savefile = [savefile,'.mat'];

if (length(tsamples) == 2)
  samples = combvec(tsamples{1},tsamples{2});
else
  samples = tsamples{1};
end

rps_size = length(samples);

if ~isequal(params.mode, 'check') && ~isequal(params.mode, 'estimator')
  settings = load(fullfile(rbmatlabtemp,testdir,'settings.mat'));
  M = settings.M;
end

all_errs = zeros(size(M,2), rps_size);
all_inds = zeros(size(M,2), rps_size);
rtimes   = zeros(size(M,2), rps_size);
rrtimes  = zeros(size(M,2), rps_size);
dtimes   = zeros(size(M,2),1);

parfor mu_ind = 1:size(M,2);
  all_errs_tmp = zeros(1,rps_size);
  all_inds_tmp = zeros(1,rps_size);
  rtimes_tmp   = zeros(1, rps_size);
  rrtimes_tmp  = zeros(1, rps_size);
  tmodel       = model;
  tempsamples  = samples;
  trp          = range_params;
  U_H          = [];

  disp(['processing parameter vector ',num2str(mu_ind),...
    '/',num2str(size(M,2))]);


  if isfield(params, 'email') && (mod(mu_ind, mail_interval) == 0) ...
      && ~isempty(params.email)
    sendmail(params.email, ['processing parameter vector ', ...
                            num2str(mu_ind),...
                            '/',num2str(size(M,2))]);
  end

  tmodel = tmodel.set_mu(tmodel, M(:,mu_ind));
  if ~isequal(params.mode, 'estimator') && ...
      ~isequal(params.mode, 'check')
    tsettings    = settings;
    [sim_data,tictoc] = load_detailed_simulation(mu_ind, ...
                                                 tsettings.savepath, tmodel);
    dtimes(mu_ind) = tictoc;
    U_H = sim_data.U;
  end


  for comb = 1:rps_size
    pf_val=tempsamples(:,comb);

    for i=1:length(pf_val)
      tmodel = update_model(tmodel, trp.plot_fields{i}, pf_val(i));
    end

    modeltemp = tmodel;
    fprintf('.');

    try
      if ~isequal(params.mode, 'check')
        [rb_sim_data,tictoc] = rb_simulation_tictoc(modeltemp, reduced_data);
        rtimes_tmp(comb) = tictoc;
      else
        rb_sim_data.Delta = -1;
      end

      if isequal(params.mode, 'estimator')
        errors = rb_sim_data.Delta;
      else
        if ~isequal(params.mode, 'check')
          tic;
          rb_sim_data = rb_reconstruction(modeltemp, detailed_data, ...
                                          rb_sim_data);
          rrtimes_tmp(comb) = toc;
          errors = tmodel.error_algorithm(U_H, rb_sim_data.U, ...
                                         detailed_data.grid, modeltemp);
        else
          errors = -1;
        end
      end
      [err,ind] = max(errors);

      all_errs_tmp(comb) = err;
      all_inds_tmp(comb) = ind;

    catch ME
      warning(['catched an error: ', ME.message]);
      disp('setting error to NaN');
    end
  end
  fprintf('.');

  all_errs(mu_ind,:) = all_errs_tmp;
  all_inds(mu_ind,:) = all_inds_tmp;
  rtimes(mu_ind,:)   = rtimes_tmp;
  rrtimes(mu_ind,:)  = rrtimes_tmp;

end;

[errs,mu_inds] = max(all_errs);
inds           = all_inds(mu_inds);

if ~isequal(params.mode, 'check')
  disp(['saving temporary workspace in ', tmpfile]);
  save(fullfile(rbmatlabtemp,tmpfile));
end

if length(tsamples) == 2
  ns1 = length(tsamples{1});
  ns2 = rps_size / ns1;
  errs = reshape(errs, ns1, ns2);
  inds = reshape(errs, ns1, ns2);
  mu_inds = reshape(errs, ns1, ns2);
  rtimes = reshape(rtimes, ns1, ns2, size(M,2));
end

ratimes = mean(sum(rtimes,3));

latextable = '';

for comb = 1:length(samples)
  for i = 1:length(tsamples)
    latextable = [latextable, '$ ', pf_descr{i},...
                              '= ', num2str(samples(i,comb))];
  end
  latextable = [latextable, '$ & $', num2str(ratimes(comb)),...
                            '$ & $', num2str(errs(comb)),...
                            '$ \\hline', sprintf('\n')];
end

dattable = sprintf('Simulation\tDimensionality\tRuntime[s]\tError\n');
dattable = [dattable, sprintf('Detailed\tH=%i\t%f\t0\n',...
                               detailed_data.grid.nelements,...
                               mean(dtimes) )];
for comb = 1:length(samples)
  dimtemp = '{';
  for i = 1:length(tsamples)
    dimtemp = [dimtemp, sprintf('%s=%i', pf_descr{i}, samples(i,comb))];
    if i~=length(tsamples)
      dimtemp = [dimtemp, ','];
    end
  end
  dattable = [dattable, sprintf(['Reduced\t', dimtemp, '\t%f\t%f\n'],...
                                ratimes(comb), errs(comb))];
end

if isequal(params.mode, 'estimator')
  bound = 1000;
else
  bound = 1;
end;

title('time index of maximum l2-error');
%set(gca,'Zscale','log');


if isfield(params, 'email') && ~isempty(params.email);
  sendmail(params.email, 'computation finished! ');
end

output.errs       = errs;
output.inds       = inds;
output.mu_inds    = mu_inds;
output.dtimes     = dtimes;
output.rtimes     = rtimes;
output.rrtimes    = rrtimes;
output.latextable = latextable;
output.dattable   = dattable;
if ~isfield(model, 'M_by_N_ratio')
  output.M_by_N_ratio = model.M_by_N_ratio;
end
output.tsamples   = tsamples;
output.samples    = samples;
output.pf_descr   = pf_descr;
output.run_name   = params.run_name;
output.stab_limit = stab_limit;
output.bound      = bound;
output.all_errs   = all_errs;
output.all_inds   = all_inds;

save(fullfile(rbmatlabtemp, savefile));

end

function model = update_model(model, plot_field, new_value)
%function model = update_model(model, plot_field, new_value)

  if ischar(plot_field)
    model.(plot_field) = new_value;
  else
    model = plot_field(model, new_value);
  end
end

function [rb_sim_data, tictoc] = rb_simulation_tictoc(varargin)
  [rb_sim_data, tictoc] = tictoc_wrapper(@rb_simulation, varargin{:});
end
