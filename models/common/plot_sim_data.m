function p = plot_sim_data(model,model_data,sim_data,plot_params)
%function p = plot_sim_data(model,model_data,sim_data,plot_params)
% function performing the plot of the simulation results as specified in model.
%
% parameters:
%   plot_params : parameter structure controlling the output of the plot. Refer
%                 to the documentation of 'model.plot_sim_data' for more
%                 details. For time dependent problems 'plot_params' are often
%                 passed to plot_element_data() and plot_sequence().
%
% required fields of model:
%   plot_sim_data : problem specific function actually performing the plot.
%
% return values:
%   p : plot handle

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

p = model.plot_sim_data(model,model_data,sim_data,plot_params);

