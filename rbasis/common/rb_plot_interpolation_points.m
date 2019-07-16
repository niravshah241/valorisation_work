function p = rb_plot_interpolation_points(detailed_data,params)
%function p = rb_plot_interpolation_points(detailed_data,params)
%
% function plotting the empirical interpolation "magic points"
% The number of points to be plotted is specified in params

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


% Bernard Haasdonk 8.6.2008

p=figure;
if ~isfield(params,'M')
  params.M = length(detailed_data.TM{1});
end;
u = zeros(detailed_data.grid.nelements,1);
u(detailed_data.TM{1}(1:params.M)) = 1:params.M;
params.clim = [0,params.M];
plot_element_data(detailed_data.grid,u,params);
title('Interpolation points/DOFS');
colorbar;

