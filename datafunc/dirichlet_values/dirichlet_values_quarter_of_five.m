function Udirichlet = dirichlet_values_quarter_of_five(glob, params)
% function Udirichlet = dirichlet_values_quarter_of_five(glob, params)
%
% function computing dirichlet-values by pointwise evaluations
%
% required fields of params:
%   dir_middle         : x coordinate where left and right regions are separated
%   c_dir_left         : left dirichlet value
%   c_dir_right        : right dirichlet value
%
% optional fields of params:
%   mu_names : names of fields to be regarded as parameters in vector mu
%
% in 'coefficients' mode, the parameters in brackets are empty

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


% glob column check
if params.debug
  if ~isempty(glob) && size(glob,1) < size(glob,2)
    warning('coordinates in variable glob are given row-wise, but expected them to be column-wise');
    if params.debug > 2
      keyboard;
    end
  end
end
decomp_mode = params.decomp_mode;


if decomp_mode == 2
  Udirichlet = [params.c_dir_left, params.c_dir_right];
else
  i = sqrt((glob(:,1).^2 + glob(:,2).^2)) <=0.1;
  j = sqrt((1-glob(:,1)).^2 + (1-glob(:,2)).^2) <= 0.1;
  if decomp_mode == 1
    Udirichlet = {i, (1-i)};
  elseif decomp_mode == 0
    Udirichlet = params.c_dir_left * i + j * params.c_dir_right;
  end
end

%| \docupdate 
