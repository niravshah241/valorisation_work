function dofs = ldg_l2project(func,qdeg,grid,params)
%function dofs = ldg_l2project(func,qdeg,grid,params)
%
% function performing an l2 projection of an analytic function func to
% the ldg space. A quadrature of degree qdeg is used.
% 'params.dimrange' and 'params.pdeg' specify the discrete function space.
%
% func is a function to be called with local coordinates 
% 'res = func(einds, locs ,grid,params)'
% 'res' is a matrix of 'size(locs,2) x dimrange'
%
% using the notation as explained in ldgdiscfunc.m, the projection computes
% df such that <df - func, `\phi_j`> = 0 for all global basis
% functions `\phi_j`
%
% This yields (using simple domain transformations, determinant
% cancles out by division)
%
% ``\mbox{dof}(j) = \frac{<\mbox{func}, \phi_j>}{|det(F_T(j)x)|} =
%                   \int_{\hat T} \mbox{func}(F(\hat x)) \hat \phi_i(j) (\hat x)``

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

% simply pass func_phi_product and its parameters to quadrature routine.

% fast implementation:
params.evaluate_basis = @ldg_evaluate_basis;
dofs = triaquadrature(qdeg,@func_phi_product,func,grid, ...
		      params);
%dofs = dofs(:);

% generic implementation:
%params.element_quadrature = @triaquadrature;
%params.local_mass_matrix = @ldg_local_mass_matrix;
%dofs = l2project(func,qdeg,grid,params);

%| \docupdate 
