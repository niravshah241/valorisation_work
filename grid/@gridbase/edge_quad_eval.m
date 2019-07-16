function F = edge_quad_eval(grid,elids,edgeids,degree,FF)
%function F = edge_quad_eval(grid,elids,edgeids,degree,FF)
% Compute an edge integral of a scalar function in various edges
% simultaneously. Approximation by Gauss-quadratures are performed.
%
% Computes quadratures approximating integrals
%  `` \int_{e(i_k,j_k)} f(s) ds ``
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

F = edge_quad_eval_mean(grid,elids,edgeids,degree,FF);
li = sub2ind(size(grid.ECX),elids(:),edgeids(:));
F = F(:).*grid.EL(li);

