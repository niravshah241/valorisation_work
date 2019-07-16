function Xdel = delzerocolumns(X,epsilon,A)
%function Xdel = delzerocolumns(X[,epsilon,A])
% function deleting zero-columns from matrix 'X'.
%
% Detection by `norm^2<\varepsilon`. Optionally, the inner-product-matrix 'A' can be
% additionally passed for correct computation of the inner product.
%
% Parameters:
%  X: matrix whose zero columns shall be eliminated
%  epsilon: tolerance value `\varepsilon`
%  A: optional inner product matrix used for norm computation (Default: `1`)
%
% Return values:
%  Xdel: matrix whose zero columns are deleted

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


% Bernard Haasdonk 20.7.2006

if nargin<3
  A = 1;
end;

if nargin<2 || isempty(epsilon)
  epsilon=eps;
end;

i = find(sum(X.*(A*X))>epsilon);
if ~isempty(i)
  Xdel = X(:,i);
else
  Xdel=zeros(size(X,1),0);
end;

