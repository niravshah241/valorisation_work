function [L_I_diff_jac, bdir_I_diff_jac] = fv_frechet_operators_diff_implicit_gradient(...
    model, model_data, U, NU_ind)
%function [L_I_diff, bdir_I_diff] = ...
%    fv_operators_diff_implicit_gradient(model,model_data, U, ...
%                                        [NU_ind])
% computes a jacobian of implicit non-linear diffusion contributions to time
% evolution matrices at a point 'U'.
%
% function computing the jacobian in 'U' of the implicit diffusion contribution
% of `L_I` and `b_I` to the time evolution matrices for a finite volume time
% step
% `L_I U^{k+1} = L_E U^k + b_E + b_I`.
% With the help of the returned Jacobian 'L_I_diff_jac' the Frechet derivative
% `DL_I ({U})` is approximated.
%
% \note The *_implicit functions perform a 'dt' increase in 'model' \em before
% evaluating the data functions.
%
% Return values:
%  L_I_diff_jac   : a sparse matrix with jacobian of diffusion contributions to
%                   `L_I`.
%  bdir_I_diff_jac: and offset vector containing the Dirichlet value
%                   contributions of the diffusion parts.
%
% Note: This only works when diffusivity is averaged on edge points by
% arithmetic mean, adapt the function, if you want to use geometric
% or harmonic mean functions.

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

tmpmodel = model;

if nargin < 4
  NU_ind = [];
end
if isempty(NU_ind)
  tmpNU = 1:length(U);
else
  tmpNU = NU_ind';
end

tmpmodel.diffusivity_ptr = model.diffusivity_derivative_ptr;

clear_all_caches;
[LU_I_diff, bdirU_I_diff] = fv_operators_diff_implicit_gradient(tmpmodel, ...
                              model_data, U, NU_ind);
clear_all_caches;
[n,m] = size(LU_I_diff);
% Note: This only works when diffusivity is averaged on edge points by
% arithmetic mean, adapt the next two lines, if you want to use geometric
% or harmonic mean functions.
tmp = LU_I_diff * U;
L_I_diff = sparse([1:n,2:n,1:(n-1)],...
                  [tmpNU(1:n), tmpNU(1:(n-1)), tmpNU(2:n)],...
                  [1/4 * tmp', 1/2 * tmp(2:n)', 1/4*tmp(1:(n-1))'],...
                  n, m);
%L_I_diff = spdiags(LU_I_diff * U,0,n,n);
bdir_I_diff = bdirU_I_diff;

tmpmodel.diffusivity_ptr = model.diffusivity_ptr;
clear_all_caches;
[LV_I_diff, bdirV_I_diff] = fv_operators_diff_implicit_gradient(tmpmodel, ...
                              model_data, U, NU_ind);
clear_all_caches;

bdir_I_diff_jac = bdir_I_diff + bdirV_I_diff;
L_I_diff_jac = L_I_diff + LV_I_diff;


