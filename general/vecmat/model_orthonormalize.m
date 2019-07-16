function onvec = model_orthonormalize(model, model_data, vec,epsilon)
% function onvec = model_orthonormalize(vec[,A,epsilon])
%
% deprecated function. simply forward to Gram-Schmidt orthonormalization

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

disp('Orthonormalize deprecated. Use orthonormalize_[qr/gram_schmidt] instead.')

if nargin < 4
  epsilon = 1e-7; % => for nonlinear evolution this value is required.
end;

onvec = model_orthonormalize_gram_schmidt(model, model_data, vec,epsilon);

