function p = plot_vertex_data(grid,data,params)
%function p = plot_vertex_data(grid,data[, params])
% plot method for vertex data on a 2D polygonal grid.
%
% For example a `P1` functions routine can be used for triangular and
% rectangular grids.  A patch plot is performed as default.  In case of
% rectangular elements, no true bilinear interpolation is performed, but the
% patch is divided into triangles.
%
% Parameters:
%  data:   a vector of length 'grid.nvertices' with nodewise scalar values.
%  params: optional structure holding fields controlling the plot output.
%
% Return values:
%  p: This is the list of handles to the graphics primitives
%
% optional fields of params:
%   shrink_factor : if this flag is given, the elements are plotted shrinked
%   axis_equal    : if this flag is set, set axis to equal scale
%   no_lines      : if this flag is set, no lines are drawn
%   show_colorbar : if this flag is set, a colorbar is drawn (default 1)
%   colorbar_location : string specifying the position of the
%                 colorbar, e.g. 'South','EastOutside' (default), etc.
%   clim          : if this 2-vector is set, the colorbar is set to
%                   these values

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

if nargin<3
  params = [];
end;

if ~(isfield(params,'shrink_factor'))
  params.shrink_factor = 1.0;
end;

if ~(isfield(params,'axis_equal'))
  params.axis_equal = 0;
end;

if ~(isfield(params,'no_lines'))
  params.no_lines = 0;
end;

if ~(isfield(params,'show_colorbar'))
  params.show_colorbar = 1;
end;

if ~(isfield(params,'colorbar_location'))
  params.colorbar_location = 'EastOutside';
end;

if (length(data)~=grid.nvertices)
  error('length of data does not match number of elements!');
end;

nneigh = grid.nneigh;

% compute vertex coordinates and scale
XX = grid.X(grid.VI(:));        
XX = reshape(XX,size(grid.VI)); % nelements*nneigh matrix
YY = grid.Y(grid.VI(:));       
YY = reshape(YY,size(grid.VI)); % nelements*nneigh matrix

CXX = repmat(grid.CX(:),1,nneigh);
CYY = repmat(grid.CY(:),1,nneigh);

% scale coordinates
XX = (XX - CXX) *params.shrink_factor + CXX;
YY = (YY - CYY) *params.shrink_factor + CYY;

%set patch colors 
CC = data(grid.VI(:));        
CC = reshape(CC,size(grid.VI)); % nelements*nneigh matrix

if isfield(params, 'plot_mode') && isequal(params.plot_mode, '3d')
  p = patch(XX',YY',CC',CC');
  if isfield(params, 'view')
    view(params.view);
  else
    view(3);
  end
else
  p = patch(XX',YY',CC');
end

if params.axis_equal
  axis equal;
  axis tight;
end;

if params.no_lines
  set(p,'linestyle','none');
end;

if params.show_colorbar
  if isfield(params,'clim')
    set(gca,'Clim',params.clim)
  end;
  colorbar(params.colorbar_location);  
end;

