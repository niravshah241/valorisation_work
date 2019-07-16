function [flux_lin, lambda] = conv_flux_forward_difference(glob, U, params)
% function [flux_lin, lambda] = conv_flux_forward_difference(glob, U, params)
% function computing the derivative of a convective flux by forward
% difference.
%
% Return values:
% 'flux_lin' : is a matrix which row entries represent the velocity vectors
%              in the edge midpoints.
% 'lambda' : is a bound such that
%             ``\lambda \cdot \sup_u n_{jl} \cdot f'(u) \leq 1``
%            e.g. `\lambda := \frac{1}{\sup|v(x,y)|}` for `f(u) = v \cdot u`.
%            This is value only reasonable in 'decomp_mode==0', otherwise an
%            empty variable is returned.
%
% Parameters:
%  U      : a vector with evaluations of a solution `u` which are passed as
%           an argument to the flux function `f`
%  glob   : a matrix of row vectors for each coordinate dimension of the grid
%           defining the coordinates where the flux function is evaluated, in
%           case it is space dependent, i.e. we have something like `f(u,x)`.
%  params : a structure with model parameters
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

epsilon = max(abs(U)) * 1e-2;
if epsilon < eps; % U completely 0
  epsilon = 1e-2;
end;

% glob column check
if params.debug
  if ~isempty(glob) && size(glob,1) < size(glob,2)
    warning('coordinates in variable glob are given row-wise, but expected them to be column-wise');
    if params.debug > 2
      keyboard;
    end
  end
end

[fluxU,lambda] = params.conv_flux_ptr(glob,U,params);
fluxU2         = params.conv_flux_ptr(glob,U+epsilon,params);

flux_lin = (fluxU2-fluxU)/epsilon;


