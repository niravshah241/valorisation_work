function Xon = orthonormalize_qr(X,A,epsilon)
% function Xon = orthonormalize_qr(X[,A,epsilon])
%
% orthonormalization of vectors in columns of 
% matrix X to Xon by qr-decomposition of vectors and 
% the inner product matrix A can be specified additionally, i.e.
% <x,x> = x' * A * x. Then the Cholesky factoriziation of A is 
% performed A = R_M' * R_M
%
% Then a QR-decomposition of R_M * X is performed which yields the
% desired new vectors

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


% Bernard Haasdonk 22.8.2008

if nargin < 3
  epsilon = 1e-7; % => for nonlinear evolution this value is required.
end;

if isempty(X)
  Xon = zeros(size(X));
  return;
end;

if nargin<2
  A = 1;
end;

% check on identity of vectors (can happen, that numerics does not detect
% this afterwards due to rounding errors!!)

for i=1:(size(X,2)-1)
  for j=(i+1):size(X,2)
    if isequal(X(:,i),X(:,j))
      X(:,j) = 0;
    end;
  end;
end;

Xon = X;

% incomplete cholesky of inner-product matrix:
R_M = cholinc(sparse(A),'inf');

% qr decomposition of R_M * X with permutation indices E
[Q, R, E]= qr(R_M * Xon,0 );

% sort such that linear independent first and Q * R = R_M * Xon
Xon = Xon(:,E);

% search nonvanishing diagonal entries of R
ind = find(abs(diag(R))>epsilon);
Xon = Xon(:,ind);
Rind = R(ind,ind); 
Xon = Xon(:,ind) / Rind;

%disp('check orthonormality!!');
%keyboard;

% eliminate zero-columns
%nsqr = sum(Xon.^2);
%i = find(nsqr > 0.1);
%Xon = Xon(:,i);

%resort zero-columns to back
%nsqr = sum(Xon.^2);
%[s,i] = sort(-nsqr);
%Xon = Xon(:,i);
%| \docupdate 
