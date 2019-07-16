function res = cubequadrature(poldeg, func, varargin)
%function res = cubequadrature(poldeg, func, varargin)
% integration of function func over reference cube == unit cube.
% by Gaussian quadrature exactly integrating polynoms of poldeg.
%
% arguments:
%  poldeg:  degree of polynomial the quadrature rule exactly approximates
%           (0-23)
%  dim:     dimension of cube on which we want to integrate (>=2)
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


[points, weights]= get_quadrature_weights_1d(poldeg);

n = length(points);

weightss = weights';
pointss  = points;

dim = 2;

for i = 1:dim-1
  m = size(pointss,1);
  weightss = kron(weightss, weights);
  weightss = reshape(weightss, m*n, 1);
  pointss  = [repmat(pointss, n, 1), reshape(repmat(points,1,m)',n*m,1)];
end

res = quadrature(weightss,pointss,func,varargin{:});

