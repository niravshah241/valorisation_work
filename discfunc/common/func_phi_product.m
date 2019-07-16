function res = func_phi_product(lcoord,func,grid,params)
%function res = func_phi_product(lcoord,func,grid,params)
%
% auxiliary function computing the product, which is used in l2
% projecting a given analytical function. 
% lcoord is a 2-vector of local coordinates, func the analytical
% function with calling syntax func(einds,loc1,loc2,grid,params)
% res is a vector of projected values on all elements. This
% function is used for integration and results in the dof values 
% of the l2-projected function.
% params.evaluate_basis is supposed to be a function doing
% point-evaluation of basis functions

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

% Phi will be a ndofs_per_element x dimrange matrix
Phi = params.evaluate_basis(lcoord,params);
%Phitrans = Phi';
nel = grid.nelements;
% F will be a nelements x dimrange matrix
%F = func(1:nel,ones(nel,1)*lcoord(:)',grid,params);
F = func(1:nel,lcoord(:)',grid,params);
%F = func(1:nel,lcoord(:)*ones(1,nel),grid,params);
% small loop over dofs
%res = zeros(grid.nelements,params.ndofs_per_element);
%res = zeros(params.ndofs_per_element,grid.nelements);

res = Phi * F';

%for i = 1:size(Phitrans,2);
%  res(:,i) = F * Phitrans(:,i);  
%end;
%for i = 1:size(Phi,2);
%  res(i,:) = (Phi(:,i)' * F);  
%end;

%| \docupdate 
