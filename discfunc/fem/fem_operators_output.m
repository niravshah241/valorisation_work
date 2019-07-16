function v = fem_operators_output(model,model_data)
%function v = fem_operators_output(model,model_data)
%
% function computing the output vectors of an elliptic problem with
% finite element discretization

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


% B. Haasdonk 22.2.2011

disp('output to be implemented, currently dummy implementation!');

% dummy output: mean over computational domain
vec = model_data.df_info.l2_inner_product_matrix * ...
      ones(model.df_info.ndofs,1);
if model.decomp_mode == 0 % == complete
  v = vec;
elseif model.decomp_mode == 1 % == components
  v ={vec};
else % == 2 coefficients
  v = 1;
end;
