function onvec = orthonormalize(vec,A,epsilon,method)
% function onvec = orthonormalize(vec[,A,epsilon,method])
%
% orthonormalization of vector a wrt 'A'-scalarproduct. Epsilon can
% be set and as methods ''gram-schmidt' (default) and 'qr' are supported.

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

if nargin < 4
  method = 'gram-schmidt';
end;

if nargin < 3 || isempty(epsilon)
  epsilon = 1e-7; % => for nonlinear evolution this value is required.
end;

if nargin<2 || isempty(A)
  A = 1;
end;

switch method
 case 'gram-schmidt'
  warning('gram-schmidt used in orthonormalize might be inaccurate!');
  onvec = orthonormalize_gram_schmidt(vec,A,epsilon);
 case 'qr'
  onvec = orthonormalize_qr(vec,A,epsilon);
 otherwise
  error('orthonormalization method unknown');
end;
%| \docupdate 
