function [L_I,b_I] = fv_operators_implicit(model, model_data, NU_ind)
%function [L_I,b_I] = fv_operators_implicit(model, model_data[ ,NU_ind])
% 
% function computing the implicit space-discretization matrix and
% offset vector for a finite volume discretization.
% Result is a sparse matrix `L_I` and an offset vector `b_I`
% Function supports affine decomposition, i.e. different operation modes
% guided by optional field decomp_mode in params. See also the
% contents.txt for general explanation.
%
% Note that this `L_I` matrix follows the notation of the
% nonlinear-FV-draft, i.e. it is NOT the time-step-operator of the 
% linear-FV-preprint in the routine fv_operators_. In the
% linear-paper the corresponding quantity would be `\bar L_I`.
% Similarly, the `b_I` has a different sign than `b` in the linear
% paper

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


% Bernard Haasdonk 16.5.2007

decomp_mode = model.decomp_mode;

if nargin <3
  if isfield(model_data, 'TM_local')
    NU_ind = model_data.TM_local{model_data.op_ind};
  else
    NU_ind = [];
  end
end

% discriminate implicit treatment in the operators
if model.fv_impl_diff_weight ~= 0
  [L_I_diff, bdir_I_diff]  = model.operators_diff_implicit(model, model_data, NU_ind);
else
  [L_I_diff, bdir_I_diff]  = fv_operators_zero(model, model_data, NU_ind);
end
if model.fv_impl_conv_weight ~= 0
  [L_I_conv, bdir_I_conv]  = model.operators_conv_implicit(model, model_data, NU_ind);
else
  [L_I_conv, bdir_I_conv]  = fv_operators_zero(model, model_data, NU_ind);
end
[L_I_neu,  bneu_I     ]  = model.operators_neumann_implicit(model, model_data, NU_ind);

% note the sign change in the following, as L_I *u + b_I is the
% discretization in nonlin_evol papar in contrast 
% to L_I u = ... + b_I    in lin_evol case

if decomp_mode == 2
  L_I = [ L_I_diff(:)  ; L_I_conv(:);  L_I_neu(:)];
  b_I = [ bdir_I_diff(:) ; bdir_I_conv(:); bneu_I(:)];
elseif decomp_mode == 0
  L_I = L_I_diff + L_I_conv + L_I_neu;
  b_I  = bdir_I_diff + bdir_I_conv + bneu_I;
elseif decomp_mode == 1
  L_I = [L_I_diff(:); L_I_conv(:); L_I_neu(:)];
  b_I = [bdir_I_diff(:); bdir_I_conv(:); bneu_I(:)];
else 
  error(['decomp_mode: ', model.decomp_mode, ' is unknown.']);
end;

