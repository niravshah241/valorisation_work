function res = exact_function_plaplace(glob, params)
%function res = exact_function_plaplace(glob, params)
% This is the first function from
% http://eqworld.ipmnet.ru/en/solutions/npde/npde1201.pdf
%
% This implements the function
% ``
% u(x,y,t) = \left(
%                 -\frac{\lambda p}{m} x
%                 - \frac{\lambda^2}{m} t
%                 + A
%           \right)^{\frac{1}{p}}
% ``
% solving the PDE
% ``
%  \partial_t u = - \partial_x \left( u^p \partial_x u \right)
% ``
% where '[x,y]' is given by 'glob', and the function parameters are read from
% 'params'. We use this function for an EOC test of our newton_model().
%
% parameters:
%  glob:      global coordinate vectors
%  params:    parameter specifying the function
%
% required fields of params:
%  t:            time instance at which the spatial solution is computed
%
% optional fields of params:
%  plaplace_A:      scalar specifying constant `A` (Default = 1.0)
%  plaplace_lambda: scalar specifying constant `\lambda` (Default = 0.8)
%  diff_p:          scalar specifying the exponent constant 'p' (Default = 0.5)
%  diff_m:          scalar specifying the constant 'm' (Default = 1.0)

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


  t = params.t;

  if ~isfield(params, 'plaplace_A')
    params.plaplace_A = 1.20;
  end

  if ~isfield(params, 'plaplace_lambda')
    params.plaplace_lambda = 0.8;
  end
  %A      = params.plaplace_A;
  %lambda = params.plaplace_lambda;

  if ~isfield(params, 'diff_p')
    params.diff_p = 0.5;
  end
  %p = params.diff_p;

  if ~isfield(params, 'diff_m')
    params.diff_m = 1.0;
  end
  %m = params.diff_m;

  X = glob(:,1);
  %Y = glob(:,2);

%  res = (-lambda * p/m * X - lambda^2*p/m*t + A).^(1/p);

  res = X + t;

