function res = ldg_evaluate_basis(lcoord,df)
%function res = ldg_evaluate_basis(lcoord,df)
%
% evaluation of ldg reference basis `\hat \phi_i, i=1,\ldots,m` in point
% `\hat x` = 'lcoord'. The argument df can either be a ldg
% discretefunction or also be a simple structure with the
% fields pdeg and dimrange.
% Res is a ndofs_per_element x dimrange (=m) array, each row being
% evaluation of one basis vector 
% lcoord is a 2-vector with the barycentric coordinates in the
% triangle. 
%
% the basis is an orthonormal basis wrt the reference unit
% simplex obtained by orthonormalizing
%    `1,x,y,x^2, xy, y^2`, etc.
% i.e. the power_vector2. 
%
% with V being the ldg_weight_matrix and W=V' we have
%
% @code
% (\hat phi_1, ... \hat phi_nbasefct)^T =  
%
%  (c_1(x),0     ,0     ,c_2(x),0     ,0     ,... c_k(x),     0,  0)  
%  (0     ,c_1(x),0     ,0     ,c_2(x),0     ,... 0     ,c_k(x),  0     )  
%  (0     ,0     ,c_1(x),0     ,     0,c_2(x),... 0     , 0    ,  c_k(x))^T  
%
% and c_i(x)= w_i * p(x) 
% @endcode
%
% with `w_i` = i-th row of W (i-th column of V) and p(x) the powervector
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


% Bernard Haasdonk 28.1.2009

switch df.pdeg
 case 0 
  res = sqrt(2) * eye(df.dimrange);
 case {1,2,3,4}
  % case 1:
  %  res = [sqrt(2) * eye(df.dimrange), ...
  %	 6 *(lcoord(1)-1/3) * eye(df.dimrange),...
  %	 6 / sqrt(3)*(2 * lcoord(2)+lcoord(1)-1) * eye(df.dimrange)];   
  V = ldg_basis_weight_matrix(df.pdeg);
  c = V' * power_vector2(lcoord,df.pdeg);
  e = eye(df.dimrange);
  nrep=[3,6,10,15];
  res = zeros(df.dimrange,nrep(df.pdeg)*df.dimrange);
  for i=1:nrep(df.pdeg)
    ra = (1+(i-1)*df.dimrange):i*df.dimrange;
    res(:,ra)=c(i)*e;
  end;
  % transposition in order to have one row per basis function 
  res = res'; 
 otherwise
  error('pdeg not yet supported!');
end;  

%| \docupdate 
