function p = lin_ds_plot_sim_data(model,model_data,sim_data,params)
%function p = lin_ds_plot_sim_data(model,model_data,sim_data[,params])
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

p1 = model.plot_sim_data_state(model,model_data,sim_data,params);
p2 = model.plot_sim_data_output(model,model_data,sim_data,params);
p = [p1(:); p2(:)];

