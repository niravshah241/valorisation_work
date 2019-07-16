function res = intervalquadrature(poldeg,func,varargin)
%function res = intervalquadrature(poldeg,func,varargin)
% integration of function func over reference interval == unit interval.
% by Gaussian quadrature exactly integrating polynoms of poldeg.
%
% arguments:
%  poldeg:  degree of polynomial the quadrature rule exactly approximates
%           (0-23)
%  func:    is a function getting a local coordinate vector (1d) and giving a
%           (vectorial or scalar) result
%  varargin: optional further arguments for function
%
% return values:
%  res:     result of quadrature
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


% Bernard Haasdonk 27.8.2009

[ points, weights ] = get_quadrature_weights_1d(poldeg);

res = quadrature(weights,points,func,varargin{:});

return

