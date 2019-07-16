function onvec = orthonormalize_gram_schmidt(vec,A,epsilon)
% function onvec = orthonormalize(vec[,A,epsilon])
%
% Gram-Schmidt orthonormalization of vectors in columns of 
% matrix vec to onvec. Almost zero vectors are deleted.
% the inner product matrix A can be specified additionally, i.e.
% <x,x> = x' * A * x;
%
% note: a threshold epsilon is involved in this routine for detecting
% zero columns. This is a quite sensitive quantity! Change with care, as
% many functions build on this!! Default is 1e-7

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

  
% Bernard Haasdonk 13.6.2002

if nargin < 3
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

if nargin<2
  A = 1;
end;

% check on identity of vectors (can happen, that numerics does not detect
% this afterwards due to rounding errors!!)

for i=1:(size(vec,2)-1)
  for j=(i+1):size(vec,2)
    if isequal(vec(:,i),vec(:,j))
      vec(:,j) = 0;
    end;
  end;
end;

onvec = vec;

for i = 1:size(vec,2);
  % orthogonalize next vector i and assume, that it is already
  % orthogonal to previous ones
  n = sqrt(onvec(:,i)' * A * onvec(:,i));
  if (n<epsilon) 
    onvec(:,i) = 0;
  else
    onvec(:,i) = onvec(:,i)/n;
  end;
  
  % orthogonalize remaining vectors wrt this one:
  
  A_mult_onvec_i = A*onvec(:,i);
  
  % create row-vector with projections on the created on-vector
  coeffs= A_mult_onvec_i' * onvec(:,i+1:end);
  
  % create matrix of scaled versions of actual orthonormalized vector
  coeffmat = onvec(:,i) * coeffs;
  
  % perform orthogonalization
  onvec(:,i+1:end) = onvec(:,i+1:end) - coeffmat;
  
end;

%keyboard;

% eliminate zero-columns
nsqr = sum(onvec.^2);
i = find(nsqr > 0.1);
onvec = onvec(:,i);



%resort zero-columns to back
%nsqr = sum(onvec.^2);
%[s,i] = sort(-nsqr);
%onvec = onvec(:,i);
%| \docupdate 
