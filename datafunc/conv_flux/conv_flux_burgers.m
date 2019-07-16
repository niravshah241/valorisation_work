function [flux, lambda] = conv_flux_burgers(glob, U, params)
% function flux = conv_flux_burgers(glob, params)
% function computing the convective flux of a Burgers problem.
%
% 'flux' is a '2xnpoints' matrix, representing the 'x'/'y'-coordinates of the
% velocity in the edge midpoints.
%
% required fields of params:
%  flux_vx   : x coordinate of flux vector
%  flux_vy   : y coordinate of flux vector
%  flux_pdeg : exponent `e` in Burgers term `f(u) = v \cdot u^e`

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

% i.e. constant velocity v = (params.flux_vx, params.flux_vy) times u^pdeg
% u is assumed to be bounded by 1, such that lambda can be
% specified in advance
Upower =  real(U(:).^params.flux_pdeg);
flux = [ params.flux_vx * ones(size(U(:))) .* Upower, ...
         params.flux_vy * ones(size(U(:))) .* Upower ];

umax = 1;
i = find(abs(U)>1.5,1);
if ~isempty(i) && params.verbose >=10
  disp(['U not bounded by 1 as assumed by flux, lambda may be ',...
        'too large.']);
end;
vmax = sqrt(params.flux_vy^2 + params.flux_vx^2);
lambda = 1/ (vmax * 2 * umax^(params.flux_pdeg-1));

