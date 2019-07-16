function F = edge_quad_eval_mean(grid,elids,edgeids,degree,FF)
%function F = edge_quad_eval_mean(grid,elids,edgeids,degree,FF)
% Compute an edge-average integral of a scalar function in various edges
% simultaneously. Approximation by Gauss-quadratures are performed.
%
% Computes quadratures approximating integrals
%  `` \frac{1}{|e(i_k,j_k)|} \int_{e(i_k,j_k)} f(s) ds ``
% over edges `e(i_k, j_k)` for `k=1,...,K`, where `i_k` are cell indices and
% `j_k` denotes local edge numbers.
%
% Parameters:
%  elids:    vector of length `K` of cell indices `i_k`.
%  edgeids:  vector of length `K` of local edge indices `j_k`.
%  degree:   scalar defining the degree of the quadrature rule.
%  FF:       matrix of size 'degree x grid.nedges' holding evaluations of a
%            function `f:\mathbb{R} \to \mathbb{R}` in quadrature points. This
%            is usually obtained by a call @code FF=f(edge_quad_points(elids,
%            edgeids, degree) @endcode
%
% Return values:
%  F: vector of length `K` holding the quadrature evaluations on edges
%     `e(i_k,  j_k)`.
%
% See also: edge_quad_point() for construction of quadrature points.

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


% Bernard Haasdonk 13.5.2007

li = sub2ind(size(grid.ECX),elids(:),edgeids(:));
nedges = length(elids);
switch degree
 case {0,1}
  % simple midpoint integration 
  % x0 = 0.5, w0 = 1;
  F = FF(:);  
 case 2
  % gauss quadrature of order 2   
  % x0 = 0.5-0.5* sqrt(1/3)    w0 = 0.5
  % x1 = 0.5+0.5* sqrt(1/3)    w1 = 0.5
  F = 0.5*(FF(1:nedges) + FF((nedges+1):end));
%  PP = [P1*x0 + P2*(1-x0), P1*x1 + P2*(1-x1)];
 case 3
  % gauss quadrature of order 3   
  % x0 = 0.5-0.5* sqrt(3/5)    w0 = 5/18
  % x1 = 0.5                   w0 = 8/18
  % x2 = 0.5+0.5* sqrt(3/5)    w0 = 5/18 
  F = (5*FF(1:nedges) + 8*FF((nedges+1):(2*nedges)) + ...
       5*FF((2*nedges+1):end))/18;
 otherwise
  error('quadrature degree not implemented!');
end;




