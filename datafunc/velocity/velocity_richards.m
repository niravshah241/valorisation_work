function [vel,lambda] = velocity_richards(glob, params)
%function [vel,lambda] = velocity_richards(glob, params)
%
% function evaluating a function in the list of global coordinates
% specified in the columns of glob. Result is a matrix of velocity
% vectors as columns of vel.
%
%
% Linear combination of components by coefficients then yields the
% complete evaluation.

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


% Martin Drohmann 23.9.2009

% glob column check
if params.debug
  if ~isempty(glob) && size(glob,1) < size(glob,2)
    warning('coordinates in variable glob are given row-wise, but expected them to be column-wise');
    if params.debug > 2
      keyboard;
    end
  end
end
X = glob(:,1);
vel = zeros(length(X),2);

%      p_mu    = spline_select(model);
%      [ breaks, coeffs, pieces, order ] = unmkpp(p_mu);
%      p_mu_d  = mkpp(breaks, coeffs(1:order-1) .* [order-1:-1:1]);

%      denom = 1 + ppval(p_mu, X);
[res1, res2] = inv_geo_trans_derivative(params,glob,{1,2,[1,1],[1,2]},...
                                                    {1,2,[2,1],[2,2]},2);
d1  = res1{3} + res2{3};
d2  = res1{4} + res2{4};
vel(:,1) = ( -params.k .* (res1{1} .* d1 + res1{2} .* d2) );
vel(:,2) = ( -params.k .* (res2{1} .* d1 + res2{2} .* d2) );

lambda = 0;
%| \docupdate 
