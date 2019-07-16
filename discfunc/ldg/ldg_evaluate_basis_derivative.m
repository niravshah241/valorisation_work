function res = ldg_evaluate_basis_derivative(lcoord,params)
%function res = ldg_evaluate_basis_derivative(lcoord,params)
%
% evaluation of ldg reference basis 
% derivatives `D(\hat \phi_i), i=1,...,m` in point
% `\hat x` = 'lcoord'. The argument params is a structure with fields 'pdeg' and 'dimrange'.
% Res is a cell array, 'res{i}' = `D \hat \phi_i`
% and `D \hat \phi_i` is a 'dimrange x 2' (==ndimworld) array 
% 'lcoord' is a 2-vector with the barycentric coordinates in the
% triangle
%
% with V being the ldg_weight_matrix and W=V' we have
%
% @code
% (\hat phi_1, ... \hat phi_nbasefct) =  
%
%  (c_1(x),0     ,0     ,c_2(x),0     ,0     ,... c_k(x),     0,  0)  
%  (0     ,c_1(x),0     ,0     ,c_2(x),0     ,... 0     ,c_k(x),  0     )  
%  (0     ,0     ,c_1(x),0     ,     0,c_2(x),... 0     , 0    ,  c_k(x))  
% @endcode
%
% and `c_i(x)= w_i * p(x)` 
%
% with `w_i` = i-th row of 'W' (i-th column of 'V') and `p(x)` the powervector
%
% with k = number of scalar base functions = dim('powervector2')
% 
% Hence the derivatives also have similar repeating structure.
% ``
% D \hat \phi_1 = \left[\begin{array}{cc} \multicolumn{2}{c}{w_1 \cdot  D p(x) } \\
%                                                              0 & 0 \\
%                                                              0 & 0 \end{array}\right]
% `` ``
% D \hat \phi_2 = \left[ \begin{array}{cc}                     0 & 0 \\
%                                         \multicolumn{2}{c}{w_1 \cdot D p(x) } \\
%                                                              0 & 0 \end{array}\right]
% ``
% etc.
%
% In the vectorial case, this computation is very storage redundant, as
% most of the `D\hat \phi_i` are zeros, only one nonzero line...
% definitely space for improvement.

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

switch params.pdeg
 case 0 % piecewise constant => zero derivative 
  res = cell(1,params.dimrange);
  for i=1:params.dimrange
    res{i}= zeros(params.dimrange,2);
  end;
 case {1,2,3,4}
  % case 1:
  %  res = [sqrt(2) * eye(df.dimrange), ...
  %	 6 *(lcoord(1)-1/3) * eye(df.dimrange),...
  %	 6 / sqrt(3)*(2 * lcoord(2)+lcoord(1)-1) * eye(df.dimrange)];   
  V = ldg_basis_weight_matrix(params.pdeg);
  nrep=[3,6,10,15];
  nbasefct = params.dimrange*nrep(params.pdeg);
  res = cell(1,nbasefct);
  Dp = power_vector2_derivative(lcoord,params.pdeg);
  WDp = V'*Dp;  
  for i=1:nrep(params.pdeg)
    for j = 1:params.dimrange
      res{(i-1)*params.dimrange+j}= zeros(params.dimrange,2);
      res{(i-1)*params.dimrange+j}(j,:) = WDp(i,:);
    end;
  end;  
 otherwise
  error('pdeg not yet supported!');
end;  
%| \docupdate 
