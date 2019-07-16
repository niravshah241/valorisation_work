function p = plot(grid,params)
%function p = plot(grid [,params])
% plot of a rectgrid via plot_polygon_grid()
%
% see help plot_polygon_grid() for further information
%
% @copydoc ::gridbase::plot_polygon_grid()

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

if (nargin <2)
  params = [];
end;

% simply forward the call
p = plot_polygon_grid(grid,params);

