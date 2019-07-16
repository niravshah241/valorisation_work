function [nmodel,nddata,plot_params] = renew_model(model, detailed_data)
% function [nmodel,nddata] = renew_model(model, detailed_data)
% change fields of old 'param' structure to new 'model' structure with
% excessive use of function pointers.
%
% This function can be used to make old demos work. Simply load the old
% 'params' and 'detailed_data' and then call 
% @code
%  [model,detailed_data,plot_params] = renewmodel(params, detailed_data);
%  demo_rb_gui(model, detailed_data, [],plot_params,'title');
% @endcode

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


if isequal(model.rb_problem_type, 'nonlin_evol')
  nmodel = nonlin_evol_model_default;
else
  nmodel = lin_evol_model_default;
end

nmodel.init_values_ptr      = str2func(['init_values_', model.name_init_values]);
model=rmfield(model,'name_init_values');

if isfield(model, 'name_dirichlet_values')
  nmodel.dirichlet_values_ptr = str2func(['dirichlet_values_', model.name_dirichlet_values]);
  model=rmfield(model,'name_dirichlet_values');
end

if isfield(model, 'name_neuman_values')
  nmodel.neumann_values_ptr   = str2func(['neumann_values_', model.name_neuman_values]);
  model=rmfield(model,'name_neuman_values');
end

nmodel.conv_flux_ptr        = str2func(['conv_flux_', model.name_flux]);
model=rmfield(model,'name_flux');

if ~isequal(model.name_diffusive_num_flux, 'none')
  nmodel.num_diff_flux_ptr    = str2func(['fv_num_diff_flux_', model.name_diffusive_num_flux]);
  nmodel.fv_expl_diff_weight = 1.0;
end
model=rmfield(model,'name_diffusive_num_flux');

if ~isequal(model.name_convective_num_flux, 'none')
  if isequal(model.name_convective_num_flux, 'enquist-osher')
    nmodel.num_conv_flux_ptr = @fv_num_conv_flux_engquist_osher;
  else
    nmodel.num_conv_flux_ptr    = str2func(['fv_num_conv_flux_', strrep(model.name_convective_num_flux,'-','_')]);
  end
  if ~model.flux_linear
    conv_lin_name = ['velocity_', model.name_convective_num_flux];
    if exist(conv_lin_name)
      nmodel.velocity_ptr = str2func(conv_lin_name);
    else
      warning(['using forward difference linearization for non-linear convective flux. Maybe, you need to implement:', conv_lin_name]);
      nmodel.velocity_ptr = @velocity_forward_difference;
    end;
  end
  nmodel.fv_expl_conv_weight = 1.0;
end
model=rmfield(model,'name_convective_num_flux');

nmodel.init_values_algorithm = str2func(model.init_values_algorithm);
model=rmfield(model,'init_values_algorithm');

nmodel.implicit_operators_algorithm = str2func(model.implicit_operators_algorithm);
model=rmfield(model,'implicit_operators_algorithm');

nmodel.inner_product_matrix_algorithm = str2func(model.inner_product_matrix_algorithm);
model=rmfield(model,'inner_product_matrix_algorithm');

model.bnd_rect_index = [];
model.flux_quad_degree = 1;
model.divclean_mode = 0;

plot_params = [];
[plot_params, model] = copy_and_rm(plot_params, model, {'yscale_uicontrols', 'xscale_uicontrols', 'xscale_gui', 'yscale_gui', 'axis_equal', 'clim', 'bind_to_model', 'no_lines', 'geometry_transformation', 'show_colorbar', 'colorbar_mode'});

nmodel     = structcpy(nmodel, model);

nmodel     = model_default(nmodel);

model_data = gen_model_data(nmodel);
nddata     = structcpy(detailed_data, model_data);

if isfield(detailed_data, 'BM')
  nddata.BM=cell(1);
  nddata.TM=cell(1);
  nddata.QM=cell(1);
  nddata.ei_info=cell(1);
  nddata.BM{1} = detailed_data.BM;
  nddata.QM{1} = detailed_data.QM;
  nddata.TM{1} = detailed_data.TM;
  nddata.ei_info{1} = detailed_data.ei_info;
  nddata.implicit_crb_index = 1;
  nddata.explicit_crb_index = 1;
end


end

function [np,op] = copy_and_rm(np,op,list)
  for i=1:length(list)
    if isfield(op, list{i})
      np.(list{i}) = op.(list{i});
      op=rmfield(op,list{i});
    end
  end
end
%| \docupdate 
