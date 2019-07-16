function [flux, lambda] = conv_flux_linear(glob,U,params)
%function [flux[, lambda]] = conv_flux_linear(glob,U,params)
% function computing the convective flux `f(u) = v u` of a convection
% problem.
%
% required fields of params:
%  t            : real time value in case of time-dependent flux
%  velocity_ptr : pointer to function providing the velocity field `v`.

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


% Bernard Haasdonk 4.4.2006

% determine affine_decomposition_mode as integer  
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

vel = params.velocity_ptr(glob,params);

if decomp_mode == 2 % i.e. 'coefficients'
  flux = vel; % copy of coefficient vector from linear flux
elseif decomp_mode == 0 % i.e. 'none'
  if ~isempty(vel) % otherwise size 0/0 => 0/1 !!!!!
    flux   = [vel(:,1).* U(:), vel(:,2).* U(:)];
    lambda = 1/max(sqrt(vel(:,1).^2 + vel(:,2).^2));
  else
    flux   = zeros(0,2);
    lambda = 0;
  end;
elseif decomp_mode == 1 % i.e. 'components'
  flux = cell(size(vel));
  Q = length(vel);
  UU = [U(:),U(:)];
  for q = 1:Q
    flux{q}= vel{q}.*UU;
  end;
end;

