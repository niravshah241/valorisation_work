function p = lin_ds_from_lin_evol_plot_sim_data(model,model_data,sim_data,params)
%function p = lin_ds_from_lin_evol_plot_sim_data(model,model_data,sim_data,params)
%
% plot of dynamical system simulation results. But not as a
% dynamical system trajectory but as a lin_evol
% model. Correspondingly, the model_data is expected to contain 
% model_data.base_model_data. The handle to the plot is returned.
%
% Return values:
%   p    plot handles

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


% Bernard Haasdonk 23.1.2010

if ~isfield(params,'no_lines');
  params.no_lines = 1;
end;
if ~isfield(params,'clim');
  params.clim = [0,1];
end;
lin_evol_from_ds_sim_data.U = sim_data.X;
lin_evol_from_ds_sim_data.y = sim_data.Y;
p = plot_sim_data(model.base_model,model_data.base_model_data,...
	      lin_evol_from_ds_sim_data,params);

