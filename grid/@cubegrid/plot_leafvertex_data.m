function p = plot_leafvertex_data(grid,data,params)
%function p = plot_leafvertex_data(grid,data,params)
% plot method for data on the leaf level of the grid.
%
% Currently only implemented for 2D. A patch plot is performed, such that the
% vertices have corresponding values as given by data, the values in between
% are interpolated by MATLAB.
%
% Return values:
%  p: this is the list of handles to the graphics primitives
%
% See also: plot_vertex_data() for optional params fields

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


% Bernard Haasdonk 23.3.2007

if nargin==1
  params = [];
  data = zeros(grid.nvertices,1);
end;

dim = grid.dimension;

if (dim~=2)
  error('plot_leafvertex_data not yet implemented for other than 2D');
end;

% get minimal information required for plotting
gridtmp = gen_plot_data(grid);

% forward call to 2d vertex plot routine
p=plot_vertex_data(gridtmp,data,params);

end

