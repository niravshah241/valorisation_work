function p = fv_plot(grid,dofs,params)
%function p = fv_plot(grid,dofs,params)
% routine plotting a single fv function of
% fv_functions.
%
% Simple forward to plot_element_data(). For sequences, simple call of
% plot_sequence() can make use of this function
%
% Parameters:
%   params: plot parameters see plot_element_data() for details.
%   dofs:   vector of degrees of freedom of finite volume diescrete function.
%
% Return values:
%   p: graphic handle to plot

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


% Bernard Haasdonk 3.9.2009

p = plot_element_data(grid,dofs,params);

