function model_data = lin_ds_gen_model_data(model)
%function model_data = lin_ds_gen_model_data(model)
%
% function computing components of matrices and inner product matrix

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



% Bernard Haasdonk 1.3.2010

if model.affinely_decomposed
  model.decomp_mode = 1;
  model_data.A_components = model.A_function_ptr(model,[]);
  model_data.B_components = model.B_function_ptr(model,[]);
  model_data.C_components = model.C_function_ptr(model,[]);
  % no D components
  %  model_data.D_components = model.D_function_ptr(model);
  model_data.x0_components = model.x0_function_ptr(model,[]);
end;
model_data.G = model.G_matrix_function_ptr(model,[]);

