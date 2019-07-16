function res = df_func_difference_sqr(lcoord,dofs,func,grid,params)
%function res = df_func_difference_sqr(lcoord,dofs,func,grid,params)
%
% auxiliary function computing the difference, which is used in l2
% error computation between a function allowing local element
% evaluation and a ldg function multiplied with integration_element 
% ( = * triangle volume) due to transformation formula
%
% lcoord is a 2-vector of local coordinates, dofs a dof vector of
% the discrete function, func the analytical
% function with calling syntax func(einds,locs,grid,params)
% res is a vector of projected values on all elements. This
% function can be used for integration and results in the l2error

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

nel = grid.nelements;
% F and DF will be a nelements x dimrange matrix
if isa(func, 'ldgdiscfunc')
  func.grid = grid;
  F = func(1:nel,lcoord,params);
else
  F = func(1:nel,lcoord, grid, params);
end
DF = params.evaluate(dofs,1:nel,lcoord,grid,params);

if params.dimrange == 1 
  res = (DF-F).^2;
else
  res = sum((DF-F).^2,2);
end;
res = res(:).*grid.A * 2;

%| \docupdate 
