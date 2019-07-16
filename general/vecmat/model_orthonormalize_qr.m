function onvec = model_orthonormalize_qr(model, model_data, vec,epsilon)
% function onvec = model_orthonormalize_qr(model, model_data, vec[,epsilon])
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


if nargin < 4
  epsilon = 1e-7; % => for nonlinear evolution this value is required.
end;
%epsilon = 1e-10; => for linear evolution this value was OK.
%epsilon = 1e-12; -> This causes error by returning non-orthogonal 
%                    vectors (e.g. two identical vectors)
%                    So take epsilon larger than 1e-12!

if isempty(vec)
  onvec = zeros(size(vec));
  return;
end;

A = model.get_inner_product_matrix(model_data);

% incomplete cholesky of inner-product matrix:
R_M = cholinc(sparse(A),'inf');

% qr decomposition of R_M * X with permutation indices E
%   [Q,R,E] = QR(B,0) produces an "economy size" decomposition in which E
%    is a permutation vector, so that B(:,E) = Q*R.

[Q, R, E]= qr(R_M * vec,0 );
%[Q, R]= qr(R_M * vec,0 );

%keyboard;

% sort such that linear independent first and Q * R = R_M * Xon
vec = vec(:,E);

% search nonvanishing diagonal entries of R
ind = find(abs(diag(R))>epsilon);
vec = vec(:,ind);
Rind = R(ind,ind); 
onvec = vec(:,ind) / Rind;

% eliminate zero-columns
%nsqr = sum(onvec.^2);
%i = find(nsqr > 0.1);
%onvec = onvec(:,i);

%resort zero-columns to back
%nsqr = sum(onvec.^2);
%[s,i] = sort(-nsqr);
%onvec = onvec(:,i);

% resort vectors such that original order is recovered most:

%keyboard;

K = model.inner_product(model,model_data,onvec,onvec);
if max(max(abs(eye(size(onvec,2))-K)))>1e-3
  K
  error(['check orthonormalization!! Non orthogonal vectors' ...
	 ' generated'])
end;
