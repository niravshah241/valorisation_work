function onvec = orthonormalize_old(vec,A)
% function onvec = orthonormalize_old(vec[,A])
%
% Gram-Schmidt orthonormalization of vectors in columns of 
% matrix vec to onvec. Almost zero vectors are set to zero.
% the inner product matrix A can be specified additionally, i.e.
% <x,x> = x' * A * x;

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

epsilon = 1e-10;
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

A_mult_onvec = zeros(size(onvec,1),0);
for i = 1:size(vec,2);
  % orthogonalize vector i
  for j=1:(i-1)
    %    onvec(:,i)= onvec(:,i)-(onvec(:,i)'*A*onvec(:,j))*onvec(:,j);
    onvec(:,i)= onvec(:,i)-(onvec(:,i)'*A_mult_onvec(:,j))*onvec(:,j);
  end; 
  % normalize vector i
  n = sqrt(onvec(:,i)' * A * onvec(:,i));
  if (n<epsilon) 
    onvec(:,i) = 0;
  else
    onvec(:,i) = onvec(:,i)/n;
  end;
  A_mult_onvec = [A_mult_onvec, A*onvec(:,i)];
end;

%resort zero-columns to back
%nsqr = sum(onvec.^2);
%[s,i] = sort(-nsqr);
%onvec = onvec(:,i);

% TO BE ADJUSTED TO NEW SYNTAX
%| \docupdate 
