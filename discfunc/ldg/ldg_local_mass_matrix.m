function M = ldg_local_mass_matrix(qdeg,params)
%function M = ldg_local_mass_matrix(qdeg,params)
%
% function computing the local mass matrix of a scalar ldg basis by
% numerical integration of degree qdeg. 
%
% Result is a numlocal_base_function x num_local_base_functions
% matrix with entries
%``\int_{\mbox{ref}_{\mbox{tria}}} \hat \phi_i \hat \phi_j.``
%
% The matrix should be unity due to orthonormalization of the
% basis, but due to numerical errors, this can be slightly different.

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

f = @(lcoord) gram_matrix(ldg_evaluate_basis(lcoord,params)');
M = triaquadrature(qdeg,f);
%| \docupdate 
