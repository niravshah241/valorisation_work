function [p0, p1] = get_edges(grid,gid)
%function [p0, p1] = get_edges(grid,gid)
% method determining the world coordinates of the edges of a single
% element with global index 'gid'.
%
% 'i'-th edge goes from 'p0(:,i)' to 'p1(:,i)'
%
% parameters:
%  gid: global index of element for which edge coordinates are requested
%
% return values:
%  p0: lowest coordinate of the edge
%  p1: highest coordinate of the edge
%
% @note vectorized version for multiple elements not yet implemented

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

  
% Bernard Haasdonk 27.3.2007
  
  dim = grid.dimension;
  
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

% get coordinates of elements vertices:
vi = grid.vertexindex(gid,:);
vertices = grid.vertex(vi,:);

p0 = vertices(li(:,1),:)';
p1 = vertices(li(:,2),:)';

end

