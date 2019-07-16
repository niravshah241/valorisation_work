function model_data = lin_ds_from_lin_evol_gen_model_data(model)
%function model_data = lin_ds_from_lin_evol_gen_model_data(model)
%
% function computing components of matrices and inner product matrix
% of a lin_ds model from a given lin_evol model

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


% Bernard Haasdonk 7.9.2009

if model.affinely_decomposed
  base_model = model.base_model;
  model_data.base_model_data = gen_model_data(model.base_model);
  base_model.decomp_mode = 1;
  % x0 components are u0 components
  model_data.x0_components = ...
      base_model.init_values_algorithm(base_model, ...
				       model_data.base_model_data);
  
  base_model.dt = base_model.T/base_model.nt;
  [L_I_comp, L_E_comp, b_comp] = ...
      base_model.operators_ptr(base_model, model_data.base_model_data);
%  disp('inspect operators for component extraction!');
%  keyboard;

  if (length(L_I_comp)~=1) || ~isequal(L_I_comp{1},speye(size(L_I_comp{1})))
    error('Conversion only reasonable for L_I == Identity!');
  end;
  
  % L_E = Delta t * A + Id
  if ~isequal(L_E_comp{1},speye(size(L_E_comp{1})))
    error('expected identity as first component of L_E')
  end;
  model_data.A_components = L_E_comp(2:end);
  
  model_data.B_components = b_comp;
  
  if ~isfield(base_model,'compute_output_functional')
    error('outputfunctional required in lin_evol!');
  else
    model_data.C_components = ...
	base_model.operators_output(base_model,model_data.base_model_data);
    for q=1:length(model_data.C_components)
      model_data.C_components{q} = model_data.C_components{q}';
    end;
  end;
  
  % no D components
  %  model_data.D_components = model.D_function_ptr(model);
end;

model_data.G = base_model.mass_matrix([],model_data.base_model_data.grid,[]);
%| \docupdate 
