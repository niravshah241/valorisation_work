function Udirichlet = dirichlet_values_weighted_boxes_coefficients(params)
%function Udirichlet = dirichlet_values_weighted_boxes_coefficients(params)
%
% Return is a vector of length Q with (parameter dependent) coefficients in U.
%
% values constant params.c_dir in box given by
%       (beta=1):  model.dir_box_xrange{1} and dir_box_yrange{1} 
%       (beta=0):  model.dir_box_xrange{2} and dir_box_yrange{2} 
%            for intermediate betas, the two boxes are weighted

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


% Bernard Haasdonk 27.8.2009

if ~ismember('beta',model.mu_names) 
  Udirichlet = model.c_dir;
else
  Udirichlet = model.c_dir * [model.beta; 1-model.beta];
end;
%| \docupdate 
