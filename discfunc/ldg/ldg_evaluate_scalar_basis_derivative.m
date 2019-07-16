function res = ldg_evaluate_scalar_basis_derivative(lcoord,df)
%function res = ldg_evaluate_scalar_basis_derivative(lcoord,df)
%
% evaluation of scalar ldg reference basis 
% derivatives `\nabla(\hat \phi_i), i=1,\ldots,m` in point
% `\hat x` = 'lcoord'. The argument 'df' can either be a ldg discfunc or
% a structure with fields pdeg and dimrange.
% Res is a nbasefunctions x dimworld array, 'Res(i,:)' = `\nabla \hat \phi_i`
% 'lcoord' is a 2-vector with the barycentric coordinates in the
% triangle
%
% with V being the ldg_weight_matrix and W=V' we have
%
% ``(\hat \phi_1, \ldots, \hat \phi_{\mbox{nbasefct}}) =  
%
%  (c_1(x),c_2(x),c_k(x))  ``
%
% and `c_i(x)= w_i  p(x) `
%
% with `w_i` = i-th row of W (i-th column of V) and p(x) the powervector
%
% with k = number of scalar base functions = dim(powervector2)
% 
% Hence the derivatives also have similar repeating structure.
%
% `\nabla \hat \phi_1 = [w_1 * D p(x)]' `
% etc.
%
% perhaps the general vectorial case (tensor-product basis) can be
% generalized from this version later...

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


% Bernard Haasdonk 28.8.2009

if (df.dimrange~=1)
  warning('caution: scalar basis evaluated, but df is vectorial!');
end;

switch df.pdeg
 case 0 % piecewise constant => zero gradient
  res = zeros(1,2); % 1 base function scalar, constant case.
 case {1,2,3,4}
  % case 1:
  %  res = [sqrt(2) * eye(df.dimrange), ...
  %	 6 *(lcoord(1)-1/3) * eye(df.dimrange),...
  %	 6 / sqrt(3)*(2 * lcoord(2)+lcoord(1)-1) * eye(df.dimrange)];   
  V = ldg_basis_weight_matrix(df.pdeg);
  nrep=[3,6,10,15];
  nbasefct = df.dimrange*nrep(df.pdeg);
  %res = zeros(2,nbasefct);
  Dp = power_vector2_derivative(lcoord,df.pdeg);
  res = Dp' * V;  
  res = res';
  %  WDp = V'*Dp;  
  %  for i=1:nrep(df.pdeg)
  %    for j = 1:df.dimrange
  %      res{(i-1)*df.dimrange+j}= zeros(df.dimrange,2);
  %      res{(i-1)*df.dimrange+j}(j,:) = WDp(i,:);
  %    end;
  %  end;  
 otherwise
  error('pdeg not yet supported!');
end;  
%| \docupdate 
