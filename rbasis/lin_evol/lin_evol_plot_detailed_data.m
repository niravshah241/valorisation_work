function lin_evol_plot_detailed_data(params, detailed_data,plot_params)
%function lin_evol_plot_detailed_data(model, detailed_data, plot_params)
%
% plot the reduced basis and generation information if available
%
% required fields of detailed_data:
%      RB : a matrix of RB DOF vectors
%
% if the field RB_info is given with generation-information, this
% is also plotted

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


% Bernard Haasdonk 23.5.2007

% plot reduced basis
plot_params.title = 'Orthonormal reduced basis';
plot_params.axis_equal = 1;
plot_params.clim_tight = 1; % color limits for each plot are set tight
plot_params.plot_function = @plot_element_data;
plot_data_sequence(params, detailed_data.grid,...
              detailed_data.RB, plot_params);

% plot generation information
if isfield(detailed_data,'RB_info')
  disp(detailed_data.RB_info)

  % mu-frequency
  if isfield(detailed_data.RB_info,'mu_sequence')
    plot_mu_frequency(detailed_data.RB_info.mu_sequence,params);
  end;

  % training error-indicator
  if isfield(detailed_data.RB_info,'max_err_sequence')
    figure,plot(detailed_data.RB_info.max_err_sequence);
    title('maximum training error indicator')
    ylabel([params.RB_error_indicator]);
    xlabel('basis vector number N');
    set(gca,'Yscale','log');
  end;

  % computation time
  if isfield(detailed_data.RB_info,'toc_value_sequence')
    figure,plot(detailed_data.RB_info.toc_value_sequence);
    title('computation time')
    ylabel('time [s]');
    xlabel('basis vector number N');
  end;

  % test assessment
  if isfield(detailed_data.RB_info,'max_test_error_sequence') && ...
        isfield(detailed_data.RB_info,'max_test_estimator_sequence')
    figure;
    estimators = detailed_data.RB_info.max_test_estimator_sequence(:);
    errors = detailed_data.RB_info.max_test_error_sequence(:);
    if isempty(isnan(estimators)) && ...
          isempty(isnan(errors))
      plot([estimators, errors]);
    else
      plot([estimators, errors],'x');
    end;
    legend({'estimator','true error'});
    title('maximum test error and indicator');
    ylabel('error and estimator');
    xlabel('basis vector number N');
    set(gca,'Yscale','log');
%    disp('please correct plotting in rb_lin_evol_plot_detailed_data!');
%    keyboard;
  end;

end;
 
