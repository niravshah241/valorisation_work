function [flux, lambda] = conv_flux_brooks_corey(glob, U, params)
% function [flux, lambda] = conv_flux_brooks_corey(glob, U, params)
% convective flux for Buckley-Leverett problem with Brooks-Corey functions
%
% function computing the nonlinear convective flux of a Buckley-Leverett with Brooks-Corey approximation
% problem.
% ``f(x,u) = \frac{\lambda_w(u)}{\lambda_w(u)+\lambda_n(u)} \quad 0\leq u \leq 1``.
%
% ``\lambda_w(u) = \frac{u^{\frac{2+3\lambda}{\lambda}}}{\mu_1},``
% ``\lambda_n(u) = \frac{(1-u)^2 (1-u^{\frac{2+\lambda}{\lambda})}}{\mu_2},``
% 
% required fields of params:
%      bl_lambda               : mobility factor `\lambda`
%      bl_mu1                  : viscosity of wetting phase `\mu_1`
%      bl_mu2                  : viscosity of non-wetting phase `\mu_2`
%                                information
%
% See also conv_flux_brooks_corey_derivative().
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

%X = glob(:,1);
%Y = glob(:,2);

ld       = params.bl_lambda;
mu1      = params.bl_mu1;
mu2      = params.bl_mu2;

lambda1  = U(:).^((2+3*ld)/ld)/mu1;
lambda2  = (1-U(:)).^2.*(1-U(:).^(2/ld+1))/mu2;

flux = lambda1 ./ (lambda1 + lambda2);


flux = [ flux, flux ] .* [0.3*ones(size(flux)),zeros(size(flux))];

lambda = 1/max(flux(:));

if params.debug && ( max(U)-eps > 1 || min(U) + eps < 0 )
  error('U is outside admissable bounds [0,1]');
end

if params.decomp_mode>0
  error('function is nonlinear and does not support affine decomposition!');
end

