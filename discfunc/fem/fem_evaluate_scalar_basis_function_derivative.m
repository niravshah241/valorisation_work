function res = fem_evaluate_scalar_basis_function_derivative(df,lcoord,i)
%function res = fem_evaluate_scalar_basis_derivative(df,lcoord,i)
%
% evaluation of `i`-th scalar fem reference basis
% derivatives `\nabla(\hat \phi_i)` in point
% `\hat x` = 'lcoord'. The argument 'df' can either be a ldg discfunc or
% a structure with fields pdeg and dimrange.
% Res is a 1 x dimworld array, 'Res(1,:)' = `\nabla \hat \phi_i`
% 'lcoord' is a 2-vector with the barycentric coordinates in the
% triangle
%
% with V being the fem_weight_matrix we have
%
% ``(\hat \phi_1, \ldots, \hat \phi_{\mbox{nbasefct}}) =  
%
%  (c_1(x),c_2(x),c_k(x))  ``
%
% and `c_i(x)= v_i  p(x) `
%
% with `v_i` = i-th row of V and p(x) the powervector
%
% with k = number of scalar base functions = dim(powervector2)
% 
% Hence the derivatives also have similar structure.
%
% `\nabla \hat \phi_1 = [w_1 * D p(x)]' `
% etc.
%
% perhaps the general vectorial case (tensor-product basis) can be
% generalized from this version later...
%
% for a simultaneuous computation of all gradients, see 
% fem_evaluate_scalar_basis_gradient

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


% Bernard Haasdonk 12.1.2011

if (df.dimrange~=1)
  warning('caution: scalar basis evaluated, but df is vectorial!');
end;

%switch df.pdeg
%case
  V = fem_basis_weight_matrix(df.pdeg);
%  v = V(ceil(i/df.dimrange),:);
  v = V(i,:);
%  nrep=(df.pdeg+1)*(df.pdeg+2)/2; %[3,6,10,15];
%  nbasefct = nrep;
  %res = zeros(2,nbasefct);
  Dp = power_vector2_derivative(lcoord,df.pdeg);
  res = v * Dp;  
%end;  
%| \docupdate 
