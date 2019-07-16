function p = lin_ds_plot_sim_data_state(model,model_data,sim_data,params)
%function p = lin_ds_plot_sim_data_state(model,model_data,sim_data[,params])
%
% plot of dynamical system simulation results
% a 3d plot of the coordinates specified in model.plot_indices
% is performed, e.g. plot_indices = [1,2,3]; performs plot of 
% first three coordinates
%
% Return values:
%  - p: plot handles

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


% Bernard Haasdonk 7.9.2009

if ~isfield(params,'plot_indices')
  plot_indices = 1:3;
else
  plot_indices = params.plot_indices;
end;

if ~isfield(params,'plot_title')
  plot_title = 'trajectory';
else
  plot_title = params.plot_title;
end;

figure, 
p = plot3(sim_data.X(plot_indices(1),:),...
	   sim_data.X(plot_indices(2),:),...
	   sim_data.X(plot_indices(3),:)...
	   ); title(plot_title);
axis equal;

