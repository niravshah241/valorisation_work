function [flux, lambda] = conv_flux_burgers_parabola(glob, U, params)
% function [flux, lambda] = conv_flux_burgers_parabola(glob, params)
% convective flux of a Burgers problem
%
% function computing the convective flux of a Burgers problem with parabolic
% velocity field times `u^2` rotated by an arbitrary angle around the point
% '(0.5,0.5)'.
%
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


P = [glob(:,1)'-0.5; glob(:,2)'-0.5];
Rinv = [cos(-params.vrot_angle), - sin(-params.vrot_angle); ...
        sin(-params.vrot_angle),   cos(-params.vrot_angle) ];
RP = Rinv*P; % rotated coordinates as columns
RP(1,:) = RP(1,:) + 0.5;
RP(2,:) = RP(2,:) + 0.5;
V = zeros(size(RP));
V(1,:) = (1- RP(2,:)) .* RP(2,:) * params.vmax * 4;
V(2,:) = 0;
RV = Rinv^(-1) * V;
flux = [ RV(1,:)' .* U(:).^params.flux_pdeg, ...
         RV(2,:)' .* U(:).^params.flux_pdeg ];

umax = 1;
i = find(abs(U)>1.5,1);
if ~isempty(i)
  disp(['U not bounded by 1 as assumed by flux, lambda may be ',...
        'too large.']);
end;

lambda = 1/(params.vmax * 2 * umax^(params.flux_pdeg-1));

