% script constructing a collateral reduced basis space for localized
% space operators.
%
% required variables that need to be set:
%   - 'model'
%   - 'model_data'
%   - 'plot_params'
%
% optional variables that can be set:
%   - 'CRBfname': filename of MAT-file where the computed 'detailed_data'
%                 structure is stored. If this variable is empty, it is
%                 set to '[model.name, infix, '_CRB_interpol.mat'], where
%                 'infix' is set to either 'model.model_type' if existent
%                 or the empty string.
%
% generated variables:
%   - 'detailed_data': structure containing the collateral reduced basis
%                      spaces and information on the empirical
%                      interpolation computations
%

disp('constructing collateral reduced basis')
%load(fullfile(rbmatlabtemp, params.step0file));

old_generation_mode = model.RB_generation_mode;
model.RB_generation_mode = 'none';
tic;
detailed_data = gen_detailed_data(model, model_data);
t = toc;
model.RB_generation_mode = old_generation_mode;
detailed_data.ei_info{1}.elapsed_time = t;
if isempty(CRBfname)
  if isfield(model, 'model_type')
    infix = model.model_type;
  else
    infix = '';
  end
  CRBfname = [model.name, infix, '_CRB_interpol.mat'];
end
save(fullfile(rbmatlabresult,CRBfname),...
     'detailed_data','model');
detailed_data.RB = zeros(detailed_data.grid.nelements,1);
plot_detailed_data(model,detailed_data,plot_params);
close(gcf-2);
