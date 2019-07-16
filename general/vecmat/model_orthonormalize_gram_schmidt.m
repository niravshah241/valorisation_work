function onvec = model_orthonormalize_gram_schmidt(model, model_data, vec,epsilon)
% function onvec = model_orthonormalize_gram_schmidt(model, model_data, vec[,epsilon])
%
% Gram-Schmidt orthonormalization of vectors in columns of
% matrix vec to onvec. Almost zero vectors are deleted.
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
% Martin Drohmann 14.5.2009
warning('gram-schmidt orthonormalization might be inaccurate!');

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
  n = sqrt(model.inner_product(model, model_data, onvec(:,i), onvec(:,i)));
  if (n<epsilon)
    onvec(:,i) = 0;
  else
    onvec(:,i) = onvec(:,i)/n;
  end;

  % orthogonalize remaining vectors wrt this one:

  A_mult_onvec_i = model.get_inner_product_matrix(model_data)*onvec(:,i);

  % create row-vector with projections on the created on-vector
  coeffs= A_mult_onvec_i' * onvec(:,i+1:end);

  % create matrix of scaled versions of actual orthonormalized vector
  coeffmat = onvec(:,i) * coeffs;

  % perform orthogonalization
  onvec(:,i+1:end) = onvec(:,i+1:end) - coeffmat;

end

% eliminate zero-columns
nsqr = sum(onvec.^2);
i = nsqr > 0.1;
onvec = onvec(:,i);

%resort zero-columns to back
%nsqr = sum(onvec.^2);
%[s,i] = sort(-nsqr);
%onvec = onvec(:,i);

K = model.inner_product(model,model_data,onvec,onvec);
if max(max(abs(eye(size(onvec,2))-K)))>1e-3
  K
  error(['check orthonormalization!! Non orthogonal vectors' ...
	 ' generated'])
end;
