function p = lin_ds_plot_detailed_data(model,detailed_data,params)
%function p = lin_ds_plot_detailed_data(model,detailed_data[,params])
%
% plot of dynamical system detailed data, i.e. the reduced basis V
% and W
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

sim_data = [];
sim_data.X = detailed_data.V;
model.plot_title = 'reduced basis V';
p1 = model.plot_sim_data_state(model,detailed_data,sim_data,params);
sim_data.X = detailed_data.W;
model.plot_title = 'biorthonormal reduced basis W';
p2 = model.plot_sim_data_state(model,detailed_data,sim_data,params);
p = [p1(:);p2(:)];
