function p = plot_polygon_grid(grid,params)
%function p = plot_polygon_grid(grid[, params])
% plot method for a 2D polygonal grid. This routine can be used for triangular
% and rectangular grids.
%
% A line plot is performed as default.
%
% Parameters:
%  params: optional structure holding fields controlling the plot output.
%
% Return values:
%  p: This is the list of handles to the graphics primitives
%
% @todo For large grids, the routine can be slow. In these cases interestingly, the
% grid plotting should be implemented with patches, as that seems to be
% faster...
%
% optional fields of params:
%   color         : RGB vector of line/patch color
%   shrink_factor : if this flag is given, the elements are plotted shrinked
%   plot_patch    : if this flag is set the plot is done
%                   by colored patches
%   axis_equal    : if this flag is set, set axis to equal scale
%

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

if nargin==1
  params = [];
end;

if ~(isfield(params,'shrink_factor'))
  params.shrink_factor = 1.0;
end;

if ~(isfield(params,'plot_patch'))
  params.plot_patch = 0;
end;

if ~(isfield(params,'axis_equal'))
  params.axis_equal = 0;
end;

if ~(isfield(params,'axis_tight'))
  params.axis_tight = 0;
end;

if ~(isfield(params,'color'))
  params.color = [0,0,1];
end;

dim = 2;
nneigh = grid.nneigh;

% compute vertex coordinates and scale
XX = grid.X(grid.VI(:));        
XX = reshape(XX,size(grid.VI)); % nelements*4 matrix
YY = grid.Y(grid.VI(:));       
YY = reshape(YY,size(grid.VI)); % nelements*4 matrix

CXX = repmat(grid.CX(:),1,nneigh);
CYY = repmat(grid.CY(:),1,nneigh);

% scale coordinates
XX = (XX - CXX) *params.shrink_factor + CXX;
YY = (YY - CYY) *params.shrink_factor + CYY;

if params.plot_patch
  p = patch(XX',YY',params.color);
else
  li = [ 1:nneigh; 1:nneigh];
  li = li(:);
  li = [li(2:end);li(1)]; % => li = [1 2 2 3 3 4 4 1];
  
  XX = XX(:,li)'; % => 2nneigh X nelements matrix
  XX = reshape(XX,2,nneigh*grid.nelements);
  YY = YY(:,li)';
  YY = reshape(YY,2,nneigh*grid.nelements);
  p = line(XX,YY,'Color',params.color);
end;

if params.axis_equal
  axis equal;
end;

if params.axis_tight
% DO NOT USE axis tight here!!! takes 1 minute with R2009a !!!!!!!
  %  axis tight;
 xmin = min(XX(:));
 xmax = max(XX(:));
 ymin = min(YY(:));
 ymax = max(YY(:));
 set(gca,'Xlim',[xmin,xmax],'Ylim',[ymin, ymax]);
end;

