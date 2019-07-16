function p = plot_leafelement_data(grid,data,params)
%function p = plot_leafelement_data(grid,data,params)
% plot method for piecewise constant data on the leaf level of the grid.
%
% Currently only implemented for 2D. A patch plot is performed.
%
% Return values:
%   p: this is the list of handles to the graphics primitives
%
% see also: plot_element_data() for general possible options in params

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
  data = zeros(grid.nelements,1);
end;

dim = grid.dimension;

if (dim~=2)
  error('plot_leafelement_data not yet implemented for other than 2D');
end;

% determine elements to be plotted
ids = find(grid.isleaf);  

if isempty(ids)
  disp('no elements on current level or grid empty')
end;

% provide new grid structure, which is acceptable by plot_element_data

elid = get(grid,'leafelements');

gridtmp = cubegrid();

gridtmp.nelements = length(elid);
gridtmp.nneigh = 4;
gridtmp.VI = grid.vertexindex(elid,[1 2 4 3]);
gridtmp.X = grid.vertex(:,1);
gridtmp.Y = grid.vertex(:,2);

% get X coordinates of vertices as vector 
XX = grid.vertex(gridtmp.VI(:),1);
% reshape to matrix
XX = reshape(XX,size(gridtmp.VI));

% get Y coordinates of vertices as vector 
YY = grid.vertex(gridtmp.VI(:),2);
% reshape to matrix
YY = reshape(YY,size(gridtmp.VI));

gridtmp.CX = mean(XX,2);
gridtmp.CY = mean(YY,2);

plot_element_data(gridtmp,data,params);

end

