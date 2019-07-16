function U = fv_l2project(func,qdeg,grid,model)
%
% function performing an l2 projection of an analytical function
% func to the fv space. A quadrature of degree qdeg is used.
% model.pdeg specify the degree of the discrete fv function.
% (currently only first order supported.)

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

disp('this is superfluous, as given by l2project');

model.evaluate_basis = @fv_evaluate_basis;
U = triaquadrature(qdeg,@func_phi_product,func,grid, ...
		      model);

