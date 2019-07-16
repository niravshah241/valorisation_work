function plot_element_data_sequence(grid,data,plot_params)
%function plot_element_data_sequence(grid,data,plot_params)
%
% plot of a sequence of element_data on the given grid (constructed
% if empty). performs simple call of plot_data_sequence.

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


% Bernard Haasdonk 9.5.2007

error('deprecated! to be replaced by discfunc//common//plot_sequence!!')

plot_params.plot_function = 'plot_element_data';
plot_data_sequence(model,grid,data,plot_params);

