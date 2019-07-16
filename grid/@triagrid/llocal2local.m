function lcoord = llocal2local(grid,faceinds,llcoord)
%function lcoord = llocal2local(grid,faceinds,llcoord)
% function performing a 1D edge-local coordinate (lcoord) to
% 2D local coordinate transformation of given faces
%
% Parameters:
%  faceinds: The face indices on which the transformation shall take
%            place. '(1:3)'
%  llcoord:  single real number between 0 and 1 defining the edge-local vertex
%            coordinate.
%
% Return values:
%  lcoord: matrix of size '2 x |faceinds|' holding the local coordinates.

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


% Bernard Haasdonk 31.8.2009

% for triagrid:
%locs3 = 1-locs(1,:)-locs(2,:);

% extend cyclically for simple index arithmetic
corner_lcoord = [0,0;1,0;0,1;0,0]';
lcoord = corner_lcoord(:,faceinds)*llcoord + ...
	 corner_lcoord(:,faceinds+1)*(1-llcoord);

