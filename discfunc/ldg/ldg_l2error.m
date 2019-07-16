function l2err = ldg_l2error(dofs,func,qdeg,grid,params)
%function l2err = ldg_l2error(dofs,func,qdeg,grid,params)
%
% function computing the l2 error of an analytic function func to
% the ldg function specified by dofs. A quadrature 
% of degree qdeg is used.
%
% func is to be called with local coordinates 
% res = func(einds, loc ,grid,params)
% where res is a matrix of length(loc) x dimrange;
%
% by this also discrete function differences can be computed by
% creating a ldgdiscfunc class object from a ldg dof vector dofs2 
%
% df2 = ldgdiscfunc(dofs2,params);
% ldg_l2error(dofs,df,...);

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


% Bernard Haasdonk 5.2.2009

% simply pass func_df_difference_sqr and its parameters to quadrature routine.

%if isa(func,'ldgdiscfunc')
%  l2errssqr = triaquadrature(qdeg,@df_df_difference_sqr,df,func,grid, ...
%			  params);
%else

params.evaluate = @ldg_evaluate;
l2errssqr = triaquadrature(qdeg,@df_func_difference_sqr,dofs,func,grid, ...
			   params);
%end;
l2err = sqrt(sum(l2errssqr));
%| \docupdate 
