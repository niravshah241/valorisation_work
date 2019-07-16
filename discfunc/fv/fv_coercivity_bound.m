function alpha = fv_coercivity_bound(model)
%function alpha = fv_coercivity_bound(model)
% 
% function determining the coercivity-alpha value of the implicit
% operator of the current simulation and current parameters.
%
% currently only for implicit diffusivity with linear dependence on
% k. In more general case, extension must be performed.
%  
% required fields of model:
% alpha_over_k : alpha is assumed to depend linearly on diffusivity k
%                value can be determined by fv_estimate_coercivity_alpha.m
%                and dividing it by k

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

  
% Bernard Haasdonk 29.8.2006

% only if diffusion is discretized implicitly  
alpha = model.alpha_over_k * model.k;
  


 