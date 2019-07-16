function model = couple_N_and_M_by_c(model, c)
% function model = couple_N_and_M_by_c(model, c)
% modifies the reduced basis size fields of 'model' by a single variable.
%
% Parameters:
%  c:  sets the reduced basis sizes to
%      `(N,M) = c (N_{\mbox{max}}, c_{MbyN} N_{\mbox{max}})`. For time-adaptive
%      schemes with several collateral reduced basis spaces, the following
%      formula is applied: `(N,M^k) = c (N_{\mbox{max}}), c_{MbyN}
%      \frac{M^k_{\mbox{max}}}{\mbox{max}_{k=0,...,K} M^k_{\mbox{max}}}
%      N_{\mbox{max}}.`
%
% Required fields of model:
%  Nmax         : maximum number of reduced basis vectors
%  Mmax         : maximum number of collateral reduced basis vectors
%  M_by_N_ratio : ratio constant `c_{MbyN}` fixing the ratio between
%                 collateral and reduced basis vector numbers
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

  model.N = max(round(c*model.Nmax),1);
  maxMmax = max(model.Mmax);
  model.M = arrayfun(@(x) min(x,max(round(c*x/maxMmax*model.M_by_N_ratio*model.Nmax),1)), model.Mmax);
end
