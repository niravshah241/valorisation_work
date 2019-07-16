function res = gram_matrix(X1,X2)
%function res = gram_matrix(X1,X2)
%
% function computing a gram matrix between the two vector sets X1
% and X2

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

if nargin < 2
  X2 = X1;
end;
dim = size(X1,1);
n1 = size(X1,2);
n2 = size(X2,2);
XX1 = repmat(reshape(X1,dim,n1,1),[1,1,n2]);
XX2 = repmat(reshape(X2,dim,1,n2),[1,n1,1]);
res = reshape(sum(XX1.*XX2,1),n1,n2);
%| \docupdate 
