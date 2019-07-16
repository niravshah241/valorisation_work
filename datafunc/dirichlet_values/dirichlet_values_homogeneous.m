function Udirichlet = dirichlet_values_homogeneous(glob, params)
% function Udirichlet = dirichlet_values_homogeneous(glob, params)
% function computing homogeneous dirichlet-values by
% pointwise evaluations
%
% required fields of params:
%   c_dir         : constant dirichlet value
%

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
  Udirichlet = params.c_dir;
elseif decomp_mode == 0
  Udirichlet = params.c_dir * ones(length(glob),1);
elseif decomp_mode == 1
  Udirichlet = { ones(length(glob),1) };
end

%| \docupdate 
