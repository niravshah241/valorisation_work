function model_data = lin_evol_gen_model_data(model)
%function model_data = lin_evol_gen_model_data(model)
%
% Function generating large model data, i.e. grid, which is not to
% be stored in the model, but required for numerics
% The reason is, that the model is passed to all kind of routines,
% hence the model is forbidden to contain high-dimensional objects.

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

model_data      = [];
model_data.grid = construct_grid(model);

% I think "inner_product_matrix_algorithm" should be replaced by mass_matrix, 
% as the latter can be built generically by local_mass_matrix
% function, e.g. for ldg functions
%model_data.W    = model.inner_product_matrix_algorithm(model, model_data);
%  the following can, e.g. be set to fv_mass_matrix
model_data.W    = model.mass_matrix([], model_data.grid,[]);

