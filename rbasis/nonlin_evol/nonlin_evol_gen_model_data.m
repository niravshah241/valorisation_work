function model_data = nonlin_evol_gen_model_data(model)
% function model_data = nonlin_evol_gen_model_data(model)
%
% method which produces all H dependent data, that is required for a 
% detailed simulations and is independent of the parameter mu.
% 
% Generated fields of model_data are:
%  grid   : the grid to be used in detailed simulations
%  W      : inner product weighting matrix

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


model_data.grid = construct_grid(model);

model_data.W    = model.inner_product_matrix_algorithm(model, model_data);

%if isequal(params.model_type, 'implicit_nonaffine_linear')
%  clear_all_caches;
%  [model_data.implicit_operator, model_data.implicit_constant ] = ...
%    fv_operators_diff_implicit(model, model_data, [], []);
%  clear_all_caches;
%end
%| \docupdate 
