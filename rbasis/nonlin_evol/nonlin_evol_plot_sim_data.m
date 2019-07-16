function p = nonlin_evol_plot_sim_data(model, model_data, sim_data, plot_params)
%function nonlin_evol_plot_sim_data(model, model_data, sim_data, plot_params)
%
% function plotting simulation results

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


% Martin Drohmann 15.12.2009

  if isfield(plot_params, 'bind_to_model') ...
      && plot_params.bind_to_model
    plot_params = copy_model_data_to_plot_params(model, ...
                                                 plot_params);
  end
  %     plot_params.plot = model.plot;
  if ~isfield(plot_params,'plot')
    plot_params.plot = model.plot;
  end;
  p = plot_sequence(sim_data.U, model_data.grid, plot_params);

