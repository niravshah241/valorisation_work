function [INC, b_I] = fv_explicit_space(model,model_data,U,NU_ind)
% function INC = fv_explicit_space(model, model_data, U, [NU_ind])
%  
% function applying an FV-space-discretization operator starting from old
% values U corresponding to the geometry given in grid producing a
% new vector of elementwise scalars NU but only on for the
% subelements with numbers given in NU_ind. If NU_ind is empty, all
% new values NU are determined, i.e. length(NU) = length(U) =
% grid.nelements
%
% By this, the operator evaluation can be performed in a localized
% way, i.e. used for empirical interpolation in rb_nonlin_evol_simulation
%
% usual timestepping can be performed afterwards by (NU = Id - deltat *
% INC). 
%
% required fields of model:
%   verbose                  :   a verbosity level
%   fv_expl_conv_weight      : factor for the explicit part of the convective
%                              discretization operator of the scheme 
%   fv_expl_diff_weight      : factor for the explicit part of the diffusion
%                              discretization operator of the scheme 
%   fv_expl_react_weight     : factor for the explicit part of the reaction term
%                              of the scheme

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


% Martin Drohmann 9.12.2007 based on fv_conv_explicit_space by Bernard
% Haasdonk

weights.diff_weight  = model.fv_expl_diff_weight;
weights.conv_weight  = model.fv_expl_conv_weight;
weights.react_weight = model.fv_expl_react_weight;

INC = fv_space_operator(model,model_data,U,NU_ind,weights);

b_I = [];

 
