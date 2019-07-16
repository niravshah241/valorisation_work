function glob = local2global(grid,einds,loc,params)
%function glob = local2global(grid,einds,loc,params)
% function performing a local to global coordinate change of
% vectors of coordinate pairs.
%
% If the three vertices of a triangle are 'v1,v2,v3', then the
% global coordinate of a single point is
% @code glob = v1 + loc(:,1).*(v2-v1) + loc(:,2).*(v3-v1); @endcode
%
% Parameters:
%  loc: matrix of size `K \times 2` holding local barycentric coordinate pairs
%       for each cell index `i_k`, `k=1,...,K`.
%  einds: vector of cell indices `i_k`, `k=1,...,K`.
%
% Return values:
%  glob: global coordinate pairs '[X, Y]' with vectors 'X' and 'Y'
%        of length `K`.
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


% Bernard Haasdonk 2.2.2009

% for triagrid:
%locs3 = 1-locs(1,:)-locs(2,:);
loc3 = 1-loc(1)-loc(2);
VI1 = grid.VI(einds,1);
VI2 = grid.VI(einds,2);
VI3 = grid.VI(einds,3);
%X = grid.X(VI1).*locs(1,:)+grid.X(VI2).*locs(2,:)+grid.X(VI3).*locs3;
%Y = grid.Y(VI1).*locs(1,:)+grid.Y(VI2).*locs(2,:)+grid.Y(VI3).*locs3;
X = grid.X(VI1)'.*loc3+grid.X(VI2)'.*loc(1)+grid.X(VI3)'.*loc(2);
Y = grid.Y(VI1)'.*loc3+grid.Y(VI2)'.*loc(1)+grid.Y(VI3)'.*loc(2);
%keyboard;
glob = [X(:),Y(:)];

