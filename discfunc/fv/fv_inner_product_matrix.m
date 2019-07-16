function W = fv_inner_product_matrix(model,model_data)
%function W = fv_inner_product_matrix(model,model_data)
%
% function computing the l2 inner product matrix for fv-functions
% on the grid.
% By this function, inner products between functions can be
% computed by their DOF vectors: K = U1' * W * U2;

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


% Bernard Haasdonk 20.7.2006

grid = model_data.grid;

% currently simple diagonal weighting of element values sufficient
W = sparse(1:grid.nelements,1:grid.nelements,grid.A);

