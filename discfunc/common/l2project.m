function U = l2project(func,qdeg,grid,params)
%function U = l2project(func,qdeg,grid,params)
%
% function performing an l2 projection of an analytical function
% func to the discrete function space space. A quadrature of degree 
% qdeg is used. params.pdeg specify the degree of the discrete fv function.
% the params.evaluate_basis must be set to the correct basis
% evaluation algorithm. params.element_quadrature is assumed to
% point to the correct element quadrature rule, i.e. triagrid, etc.
%
% result is U the column vector of dofs

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


% Bernard Haasdonk 3.9.2009

% result is a grid.nelements x params.ndofs_per_element matrix

U = params.element_quadrature(qdeg,@func_phi_product,func,grid, ...
		      params);

% multiply locally with inverse reference mass matrix
%Utmp = reshape(U,params.ndofs_per_element,grid.nelements);
M = params.local_mass_matrix(qdeg,params);
U = U * inv(M)';
U = U(:);
%| \docupdate 
