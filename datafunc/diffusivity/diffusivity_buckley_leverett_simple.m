function diffusivity = diffusivity_buckley_leverett_simple(glob, U, params)
% function diffusivity = diffusivity_buckley_leverett_simple(glob, U, params)
%
% function computing an nonlinear diffusivity for the buckley leverett problem
% ``k(u) = K p_c'(u) \cdot \lambda_2(u) f(u) ``
% with
% ``p_c(u) = u^{-\lambda},``
% ``p_c'(u) = -\lambda u^{-\lambda-1},``
% ``f(u) = \frac{\lambda_w(u)}{\lambda_w(u)+\lambda_n(u)},``
% ``\lambda_w(u) = \frac{u^3}{\mu_1},``
% ``\lambda_n(u) = \frac{(1-u)^3}{\mu_2}.``
%
% parameters:
%      U                      : solution at points in 'glob'
%
% required fields of params:
%      diff_K                 : diffussion factor `K`
%      bl_mu_1                : viscosity factor `\mu_1`
%      bl_mu_2                : viscosity factor `\mu_2`
%      bl_lambda              : mobility parameter `\lambda`
%
% See also diffusivity_buckley_leverett_simple_derivative().
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
if params.newton_regularisation && ( ( max(U) > 1 ) ...% && max(U)-1e7*eps < 1 )...
                               || ( min(U) < 0 ) )%&& min(U)+1e7*eps > 0 ) )
  U(U>1) = 1;
  U(U<0) = 0;
end
  
if( params.debug && any(U<0))
  error('argument function is negative, which is disallowed here');
end

ld  = params.bl_lambda;
mu1 = params.bl_mu1;
mu2 = params.bl_mu2;

lambda1 = U.^3/mu1;
lambda2 = (1-U).^3/mu2;

mobility = lambda1 + lambda2;

%deg1 = (1-U).^(-0.5);
%deg1(isinf(deg1)) = realmax;

diffusivity.K = + params.diff_K * ld/mu1 ...
                  .* lambda2 ...
                  .* (real(U.^(2-ld)) ./ mobility);
 

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

