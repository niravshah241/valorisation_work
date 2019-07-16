function model_data = lin_stat_gen_model_data(model)
%function model_data = lin_stat_gen_model_data(model)
%
% generation of model data of a lin-stat model
% generated fields of model_data: grid and df_info

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


% Bernard Haasdonk

model_data = [];
model_data.grid = construct_grid(model);
model_data.df_info = feminfo(model,model_data.grid);
%model_data.inner_product_matrix = ...
%    model_data.df_info.h10_inner_product_matrix
