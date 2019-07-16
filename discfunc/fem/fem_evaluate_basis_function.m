function res = fem_evaluate_basis_function(df,lcoord,i)
%function res = fem_evaluate_basis_function(df,lcoord,i)
%
% evaluation of single fem reference basis `\hat \phi_i, i=1,\ldots,m` in point
% `\hat x` = 'lcoord'. 
%
% The argument df is a fem
% discretefunction 
% Res is a 1 x dimrange (=m) array.
% lcoord is a 2-vector with the local coordinates in the
% reference triangle. 
%
% the basis is a lagrange basis with respect to the lagrange-points
% counted from v1 to v2 then row-wise upward to v3
%
% the basis function is characterizes as linear combination of the
% standard monomial basis:
%    `1,x,y,x^2, xy, y^2`, etc.
% i.e. the power_vector2. 
%
% with V being the fem_weight_matrix we have
%
% @code
% (\hat phi_1, ... \hat phi_nbasefct)^T =  
%
%  (c_1(x),0     ,0     ,c_2(x),0     ,0     ,... c_k(x),     0,  0)  
%  (0     ,c_1(x),0     ,0     ,c_2(x),0     ,... 0     ,c_k(x),  0     )  
%  (0     ,0     ,c_1(x),0     ,     0,c_2(x),... 0     , 0    ,  c_k(x))^T  
%
% and c_i(x)= v_i * p(x) 
% @endcode
%
% with `v_i` = i-th row of V and p(x) the powervector
%
% with k = number of scalar base functions = dim(powervector2)
%
% For simultaneous evaluation of all basis vectors see fem_evaluate_basis

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

switch df.pdeg
 case 0 
  error('fem basis not reasonable for pdeg 0')
 otherwise
  V = fem_basis_weight_matrix(df.pdeg);
  v = V(ceil(i/df.dimrange),:); % ith row
  c = v * power_vector2(lcoord,df.pdeg);
  ncomp = mod(i-1,df.dimrange)+1;
  res = zeros(1,df.dimrange);
  res(ncomp) = c;
end;  

%| \docupdate 
