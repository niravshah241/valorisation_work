function [U,s] = svds_batch(X,k_per_part,n_parts)
%function [U,s] = svds_batch(X,k_per_part,n_parts)
%
% function computing k_per_part * n_parts POD modes of
% the matrix X and singular values in vector s
% function can be used if svds(X,k_per_part*n_parts) gives 
% memory problem. 
%
% exact orthogonality is not obtained and some zero vectors can be returned!

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


% Bernard Haasdonk 29.5.2010

k = k_per_part*n_parts;
U = zeros(size(X,1),k);
s = zeros(k,1);
n = size(X,2);

if (k>n)
  error('too many modes requested.');
end;

part_start_index = floor((0:n_parts)*(k/n_parts))+1;
found_zero_vector = 0;
for p = 1:n_parts
  if ~found_zero_vector
    ind = part_start_index(p):(part_start_index(p+1)-1);
    disp(['computing mode ',num2str(ind(end)),'...']);
    %  keyboard;
    [Unew,S] = svds(X,k_per_part);
    ind2 = ind(1:size(Unew,2));
    U(:,ind2) = Unew;
    s(ind2) = diag(S);
    disp(['orthogonalization ...']);
    % orthonormalization:
    X = X - U(:,1:ind(end))*(U(:,1:ind(end))'*X);
    if ~isequal(ind,ind2)
      found_zero_vector = 1;
    end;
  end;    
end;
%disp(['orthonormaization with qr ...']);
%U = orthonormalize_qr(U);
disp(['svds_batch finished ']);


