function res = fem_evaluate_basis(df,lcoord)
%function res = fem_evaluate_basis(df,lcoord)
%
% evaluation of fem reference basis `\hat \phi_i, i=1,\ldots,m` in point
% `\hat x` = 'lcoord'. The argument df is a fem
% discretefunction 
% Res is a ndofs_per_element x dimrange (=m) array, each row being
% evaluation of one basis vector 
% lcoord is a 2-vector with the barycentric coordinates in the
% triangle. 
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
  c = V * power_vector2(lcoord,df.pdeg);
  e = eye(df.dimrange);
  %nrep=[3,6,10,15,21,28];
  nrep = (df.pdeg+1) * (df.pdeg+2)/2;
  res = zeros(df.dimrange,nrep * df.dimrange);
  for i=1:nrep
    ra = (1+(i-1)*df.dimrange):i*df.dimrange;
    res(:,ra)=c(i)*e;
  end;
  % transposition in order to have one row per basis function 
  res = res'; 
end;  

%| \docupdate 
