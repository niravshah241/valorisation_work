function p = plot(grid,params)
%function p = plot(grid,params)
% plot method for a cubegrid.
%
% A line plot is performed. Currently implementation for more than 3D not
% implemented, as this would be a line-mess
%
% Return values:
%  p: this is the list of handles to the graphics primitives
%
% optional fields of params:
%   color : RGB vector
%   shrink_factor: if this flag is given, the elements are plotted shrinked
%   plot_level:   if this flag is nonzero, a level plot is performed
%   level  : this integer must be specified in case of level-plot    
%   plot_patch    : if this flag is set (only 2D and 3D) the plot is done
%                   by colored patches
%   axis_equal    : if this flag is set, set axis to equal scale
%
% Bernard Haasdonk 1.3.2007

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


if nargin==1
  params = [];
end;

if ~(isfield(params,'shrink_factor'))
  params.shrink_factor = 1.0;
end;

if ~(isfield(params,'plot_level'))
  params.plot_level = 0;
end;

if ~(isfield(params,'level'))
  params.level = 0;
end;

if ~(isfield(params,'plot_patch'))
  params.plot_patch = 0;
end;

if ~(isfield(params,'axis_equal'))
  params.axis_equal = 0;
end;

if ~(isfield(params,'color'))
  params.color = [0,0,1];
end;

if ~(isfield(params,'LineStyle'))
    params.LineStyle = '-'
end;

dim = grid.dimension;

% determine elements to be plotted
if (params.plot_level)
  ids = find(grid.level==params.level);
else
  ids = find(grid.isleaf);  
end;

if isempty(ids)
  disp('no elements on current level or grid empty')
end;

if ~params.plot_patch % i.e. line mode required
  % generate list of local coordinates to be connected
  % idea: plot n-1 d bottom cube, its shifted top version and the 
  % connecting lines:
  % dim = 1 =>    [1 2]
  % dim = 2 =>    [1 2; 3,4; 1 3; 2 4]
  % dim = 3 =>    [1 2; 3,4; 1 3; 2 4, 5 6; 7,8; 5 7; 6 8; 1 5 ; 2 6 ;
  %                3 7 ; 4 8]
  li = [1 2];
  for i = 2:dim
    % assume li is local coordinate li of lower n-1 d patch
    li = [li; li + 2^(i-1)];
    % now li is local coordinate li of lower and upper n-1 d patch
    li = [li; (1:2^(i-1))' ,(1:2^(i-1))'+2^(i-1)];
  end;
else % generate patch-list in local coordinates
  switch dim
   case 2
    li = [1 2 4 3];
   case 3
    li = [1 2 4 3; ...
          5 6 8 7; ...
	  1 2 6 5; ...
	  2 4 8 6; ...
	  4 3 7 8; ...
	  3 1 5 7 ];
   otherwise
    error('patch mode requires dim= 2 or 3!');
  end;
end;

Xtotal = [];
Ytotal = [];
Ztotal = [];

X = zeros(size(li,2),size(li,1));
Y = zeros(size(li,2),size(li,1));
if dim >= 3
  Z = zeros(size(li,2),size(li,1));
else
  Z = [];
end;

% stupid plotting: loop over all elements. should be vectorized sometime 

for nel = ids(:)'
  % collect globalcoordinates of points for element
  ve = grid.vertexindex(nel,:);
  co = grid.vertex(ve,:); % => co is nlocalpoints x dimension matrix
			  % midpoint of element as matrix
  cog = repmat(mean(co),size(co,1),1); 
  % scale coordinates
  co = (co - cog)*params.shrink_factor + cog;
  
  % generate coordinate lists to be plotted
  % X,Y,Z: 2 x nlines matrix
  switch grid.dimension
   case 1
    % only x-coord varying
    X(1,:) = co(li(:,1),1);
    X(2,:) = co(li(:,2),1);    
   case 2
    % x and y-coord varying
    X = reshape(co(li',1),size(X));
    Y = reshape(co(li',2),size(Y));
   case 3
    % x and y- and z-coord varying
    X = reshape(co(li',1),size(X));
    Y = reshape(co(li',2),size(Y));
    Z = reshape(co(li',3),size(Z));
  end;      
  Xtotal = [Xtotal, X];
  Ytotal = [Ytotal, Y];
  Ztotal = [Ztotal, Z];  
end;

switch grid.dimension
 case 1
  p = line(Xtotal,Ytotal,'Color',params.color,'LineStyle',params.LineStyle);  
 case 2
  if params.plot_patch
    p = patch(Xtotal,Ytotal,params.color);
  else
    p = line(Xtotal,Ytotal,'Color',params.color,'LineStyle',params.LineStyle);
  end;
 case 3
  if params.plot_patch
    p = patch(Xtotal,Ytotal,Ztotal,params.color);
  else
    p = line(Xtotal,Ytotal,Ztotal,'Color',params.color,LineStyle,params.LineStyle);
  end;
  view(3);
 otherwise
  disp('plot for dimension>3 not yet implemented!')
  return;
end;

if params.axis_equal
  axis equal;
end;

end

