function diffusivity = diffusivity_exponential(glob, U, params)
% function diffusivity = diffusivity_exponential(glob, U, params)
%
% function computing an exponential nonlinear diffusivity
% ``k(x,u) = k_0 + m u^p ``.
%
% parameters:
%      U                      : solution at points in 'glob'
%
% required fields of params:
%      diff_k0                : constant part `k_0 \in [0.1,0.5]`
%      diff_m                 : steepness factor of curve `m \in [0, 0.5]`.
%      diff_p                 : exponent `p \in [1/100, 100]`
%
% See also diffusivity_exponential_derivative().
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

diffusivity.K = params.diff_k0;
if params.diff_m ~= 0
  diffusivity.K = diffusivity.K + params.diff_m * real(real(U).^params.diff_p);
end

diffusivity.epsilon = max(diffusivity.K);

if params.debug
  if ( params.diff_m ~= 0 && params.diff_p < 1 && min(U) - eps < 0 )
                                    %(max(U)+eps > 1 || min(U) - eps < 0 ) )
    error('U is outside admissable bounds (>0)');
  end
  if ( max(abs(imag(U))) > 1e-6 )
    error('U has non-trivial imaginary addend.')
  end
end

if params.decomp_mode>0
  error('function is nonlinear and does not support affine decomposition!');
end

