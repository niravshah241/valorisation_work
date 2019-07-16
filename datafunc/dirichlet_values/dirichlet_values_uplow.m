function Udirichlet = dirichlet_values_uplow(glob, params)
% function Udirichlet = dirichlet_values_uplow(glob, params)
%
% function computing dirichlet-values by pointwise evaluations  
%
% required fields of params:
%   dir_middle         : y coordinate where upper and lower regions are separated
%   c_dir_up, c_dir_low: upper and lower dirichlet value
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
  Udirichlet = [params.c_dir_low, params.c_dir_up];
else
  i = (glob(:,2) <= params.dir_middle) ...
       + (1-(glob(:,2)-params.dir_middle)./0.2) ...
       .* (glob(:,2) >= params.dir_middle & glob(:,2) <= params.dir_middle + 0.2);
  if decomp_mode == 1
    Udirichlet = {i, (1-i)};
  elseif decomp_mode == 0
    Udirichlet = params.c_dir_low * i + (1-i) * params.c_dir_up;
  end;
end;

%| \docupdate 
