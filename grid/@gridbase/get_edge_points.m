function [P1, P2] = get_edge_points(grid,elids,edgeids)
%function [P1, P2] = get_edge_points(grid,elids,edgeids)
% function extracting edge coordinates from the grid.
%
% Given a grid of dimension `d`, this function returns the two extremal
% coordinates of edges `\{e(i_k, j_k)\}_{k=1}^K` assuming that the edge is
% `(d-1)`-dimensional cube.  Here, `i_k` are cell indices and `j_k` denotes
% local edge numbers for all `k=1,...,K`
%
% Parameters:
%  elids:    vector of length `K` of cell indices `i_k`.
%  edgeids:  vector of length `K` of local edge indices `j_k`.
%
% Return values:
%  P1:    matrix of size `K x d` holding the first  extremal coordinate of the
%         edges `e(i_k, j_k)`.
%  P2:    matrix of size `K x d` holding the second extremal coordinate of the
%         edges `e(i_k, j_k)`.

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


% Bernard Haasdonk 12.5.2007

li = sub2ind(size(grid.ECX),elids(:),edgeids(:));
edgeids_plus_1 = edgeids + 1;
i = edgeids_plus_1 > grid.nneigh;
edgeids_plus_1(i) = 1;
li_plus_1 = sub2ind(size(grid.ECX),elids(:),edgeids_plus_1(:));
P1 = [grid.X(grid.VI(li)),...
      grid.Y(grid.VI(li))];
P2 = [grid.X(grid.VI(li_plus_1)),...
      grid.Y(grid.VI(li_plus_1))];

