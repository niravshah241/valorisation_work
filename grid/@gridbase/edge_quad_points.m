function PP = edge_quad_points(grid,elids,edgeids,degree)
%function PP = edge_quad_points(grid,elids,edgeids, degree)
% get the evaluation points for a quadrature over edges of the
% given grid.
%
% Given a list of edges `\{e(i_k, j_k)\}_{i=1}^K`, where `i_k` are cell indices
% and `j_k` denotes local edge numbers, this function computes quadrature
% points
% ``\left( x_k^1, \ldots, x_k^Q \right)``
% for each `k=1,\ldots, K`.
%
% Parameters:
%  elids:    vector of length `K` of cell indices `i_k`.
%  edgeids:  vector of length `K` of local edge indices `j_k`.
%  degree:   scalar defining the degree of the quadrature rule.
%
% Return values:
%  PP:   matrix of size `K \times Q` holding the quadrature for each edge.
%        This can be used as last argument for either edge_quad_eval() or
%        edge_quad_eval_mean().
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


% Bernard Haasdonk 12.5.2007

switch degree
 case {0,1} 
  % simple midpoint integration 
  % x0 = 0.5, w0 = 1;
  li = sub2ind(size(grid.ECX),elids(:),edgeids(:));
  PP = [grid.ECX(li), grid.ECY(li)];
 case 2
  % gauss quadrature of order 2   
  % x0 = 0.5-0.5* sqrt(1/3)    w0 = 0.5
  % x1 = 0.5+0.5* sqrt(1/3)    w1 = 0.5
  x0 = 0.5 - 0.5*sqrt(1/3);
  x1 = 0.5 + 0.5*sqrt(1/3);
  [P1,P2] = get_edge_points(grid,elids,edgeids);
  PP = [P1*x0 + P2*(1-x0); P1*x1 + P2*(1-x1)];
 case 3
  % gauss quadrature of order 3   
  % x0 = 0.5-0.5* sqrt(3/5)    w0 = 5/18
  % x1 = 0.5                   w0 = 8/18
  % x0 = 0.5+0.5* sqrt(3/5)    w0 = 5/18 
  a = 0.5*sqrt(0.6);
  x0 = 0.5 - a;
  x1 = 0.5;
  x2 = 0.5 + a;
  [P1,P2] = get_edge_points(grid,elids,edgeids);
  PP = [P1*x0 + P2*(1-x0); P1*x1 + P2*(1-x1); P1*x2 + P2*(1-x2) ];
 otherwise
  error('quadrature degree not implemented!');
end;


