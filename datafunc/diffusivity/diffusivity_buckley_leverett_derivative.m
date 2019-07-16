function diffusivity = diffusivity_buckley_leverett_derivative(glob, U, params)
% function diffusivity = diffusivity_buckley_leverett_derivative(glob, U, params)
%
% function computing derivative of an nonlinear diffusivity for the Buckley-Leverett problem with
% respect to `u` : 
% `` \partial_u k(x,u) = K p_c''(u) \lambda_n(u) f_w(u) + p_c'(u) \left (\lambda_n'(u)f_w(u) + \lambda_n(u) f_w'(u) \right) ``
% `` = -K\frac{\lambda}{\mu_1} \left[
%          u^{\frac{2+2\lambda-\lambda^2}{\lambda}}
%          \left( \frac{\lambda_n(u) m'(u)}{(m(u))^2} + \frac{\lambda_n'(u)}{m(u)} \right)
%         + u^{\frac{2+\lambda-\lambda^2}{\lambda}}
%           \frac{\lambda_n(u)}{m(u)} \left(\frac{2\lambda-\lambda^2+2}{\lambda}\right) \right]``
% with
% ``p_c(u) = u^{-\lambda},``
% ``p_c'(u) = -\lambda u^{-\lambda-1},``
% ``p_c''(u) = (\lambda^2+\lambda) {u^{-\lambda-2}},``
% ``f_w(u) = \frac{\lambda_w(u)}{\lambda_w(u)+\lambda_n(u)},``
% ``\lambda_w(u) = \frac{u^{\frac{2+3\lambda}{\lambda}}}{\mu_1},``
% ``\lambda_n(u) = \frac{(1-u)^2 (1-u^{\frac{2+\lambda}{\lambda})}}{\mu_2},``
% ``\lambda_w'(u) = \frac{2+3\lambda}{\lambda} \frac{u^{\frac{2+2\lambda}{\lambda}}}{\mu_1},``
% ``\lambda_n'(u) = \frac{-1}{\mu_2} \left(2(1-u) (1-u^{\frac{2+\lambda}{\lambda})} + \frac{2+\lambda}{\lambda}(1-u)^2 (1-u^{\frac{2}{\lambda}}) \right),``
% ``m(u) = \lambda_w(u) + \lambda_n(u),``
% ``m'(u) = \lambda_w'(u) + \lambda_n'(u).``
%
% parameters:
%      U                      : solution at points in 'glob'
%
% required fields of params:
%      diff_K                 : diffussion factor `K`
%      bl_mu_1                : viscosity factor `\mu_1`
%      bl_mu_2                : viscosity factor `\mu_2`
%      bl_lambda              : material parameter `\lambda`
%
% See also diffusivity_buckley_leverett().
%
% generated fields of diffusivity:
%     epsilon: upper bound on diffusivity value
%     K: vector with diffusivity values

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

ld       = params.bl_lambda;
mu1      = params.bl_mu1;
mu2      = params.bl_mu2;

lambda1  = U.^((2+3*ld)/ld)/mu1;
lambda2  = (1-U).^2.*(1-U.^(2/ld+1))/mu2;

mobility = lambda1 + lambda2;

lambda1d = U.^((2+2*ld)/ld)/mu1*(2/ld+3);
lambda2d = -1/mu2*(2*(1-U).*(1-U.^(2/ld+1)) + (1-U).^2.*(1-U.^(2/ld)*(2/ld+1)));

mobilityd = lambda1d + lambda2d;
%gud     = 2*U./params.diff_mu1;
%hu      = lambda1+lambda2;
%hud     = gud + 2*(U-1)./params.diff_mu2;
%fu      = lambda1./(hu);
%fud     = gud./hu - (gu.*hud)./(hu.^2);

diffusivity.K = + params.diff_K * ld/mu1 ...
                 * U.^(2-ld+2/ld) .* ( ...
                 U .* ( (lambda2.*mobilityd)./mobility.^2 ...
                        + lambda2d./mobility ) ...
                 + lambda2./mobility * (2-ld+2/ld) );

diffusivity.epsilon = max(diffusivity.K);

if params.debug
  if ( any(isnan(diffusivity.K) | isinf(diffusivity.K)) )
    error('diffusion factor has invalid elements');
  elseif ( max(abs(imag(U))) > 1e-6 )
    error('U has non-trivial imaginary addend.');
  end
end

if params.decomp_mode>0
  error('function is nonlinear and does not support affine decomposition!');
end

