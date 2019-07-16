function W = inner_product_matrix(grid,params)
%function W = inner_product_matrix(grid,params)
%
% function computing the stiffness-matrix of the discrete functions
% on grid specified in params.
%
% required fields of params:
%
%     'inner_product_matrix_algorithm' : name of function, which
%         gives the matrix W for L2-norm computation U1' * W * U2
%         e.g. fv_inner_product_matrix. Arguments of this function
%         are the grid and params.

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

% currently simple diagonal weighting of element values sufficient
W = feval(params.inner_product_matrix_algorithm,grid,params);

% TO BE ADJUSTED TO NEW SYNTAX
%| \docupdate 
