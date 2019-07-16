function display(grid)
%function display(grid)
%
% default implementation for display method for unstructured polygonal grids.
% is used by rectgrid and triagrid

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
disp(['gridtype : ',class(grid)]);
disp(['number of elements : ', num2str(grid.nelements)]);
disp(['number of vertices : ', num2str(grid.nvertices)]);

%display(grid.vertex);
%display(grid.vertexindex);
%display(grid.level);
%display(grid.isleaf);

