function res = ldg_evaluate(dofs, eindices, lcoord, grid, params)
%function res = ldg_evaluate(dofs, eindices, lcoord, grid, params)
%
% method evaluating a ldg function in local coordinates in the point
% `\hat x` = lcoord in several elements simultaneously given by 
% indices eindices.
% res is a length(eindices) x dimrange vector 

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


% Bernard Haasdonk 28.1.2009

% evaluate all reference basis functions in the lcoords, 
% is a ndofs_per_element x dimrange matrix
basis_values = ldg_evaluate_basis(lcoord,params);

%res = zeros(params.dimrange,length(eindices));
%res = zeros(length(eindices),params.dimrange);
% linear combination with DOFS to get function values

d = reshape(dofs,params.ndofs_per_element,grid.nelements);

% small loop over basis functions
%for i = 1:params.ndofs_per_element
%%  res = res + basis_values(:,i) * d(i,eindices);
%  res = res + d(eindices,i)* basis_values(i,:);
%end;

%res = basis_values' * d(:,eindices); 
res = d(:,eindices)' * basis_values; 

%ldg_evaluate_basis([0,0],params)'*dofs(1,:)'
%| \docupdate 
