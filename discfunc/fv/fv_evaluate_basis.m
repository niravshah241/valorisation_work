function res = fv_evaluate_basis(lcoord,df)
%function res = fv_evaluate_basis(lcoord,df)
%
% evaluation of fv reference basis ` \hat \phi_i, i=1,...,m` in point
% `\hat x = \mbox{lcoord}`.
%
% currently only pdeg0 fv functions are supported, hence
% all values are 1.0
% Res is a dimrange x ndofs_per_element (=m) array 
% lcoord is a 2-vector with the barycentric coordinates in the
% triangle. 
%

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
  res = 1;
 case {1,2,3,4}
  error('pdeg not yet supported!');
end;  

