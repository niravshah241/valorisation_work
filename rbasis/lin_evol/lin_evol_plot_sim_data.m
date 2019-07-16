function p = lin_evol_plot_sim_data(model,model_data,sim_data,params)
%function p = lin_evol_plot_sim_data(model,model_data,sim_data.params)
%
% function plotting simulation results
%
% Optional fileds of params:
%   plot_title: Title of the plot
%   plot_output_title: title of the output plot
%   axis_equal: axis equal

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


% Bernard Haasdonk 4.9.2009

if ~isfield(params,'plot_title')
  plot_title = 'trajectory';
else
  plot_title = params.plot_title;
end;

if ~isfield(params,'plot_output_title')
  plot_output_title = 'output';
else
  plot_output_title = params.plot_output_title;
end;

% plot sequence
params.title = plot_title;
if ~isfield(params,'plot')
  params.plot = model.plot;
end
params.axis_equal = 1;
p = plot_sequence(sim_data.U,model_data.grid,params);

% plot output functional
if model.compute_output_functional
  p2 = figure;
  plot(sim_data.y);
  title(plot_output_title);
  p = [p,p2];
end;

