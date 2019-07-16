function M = fv_local_mass_matrix_rect(qdeg,model)
%function M = fv_local_mass_matrix_rect(qdeg,model)
%
% function computing the local mass matrix of a scalar fv basis by
% numerical integration of degree qdeg. 
%
% Result is a numlocal_base_function x num_local_base_functions
% matrix with entries 
% `\int_{\mbox{ref}_{\mbox{tria}}} \hat phi_i \hat phi_j`.
%
% as only first order fv elements are supported currently, 
% the "matrix" is a simply 1/2 for triangles and 1 for rectangles

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

M = 1;

