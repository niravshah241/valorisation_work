% script constructing a reduced basis space
%
% required variables that need to be set:
%   - 'CRBfname'
%   - 'plot_params'
%
% optional variables that can be set:
%   - 'detailedfname': filename of MAT-file where the computed
%                      'detailed_data' structure is stored. If this
%                      variable is empty, it is set to
%                      '[model.name, infix, '_detailed_interpol.mat'],
%                      where 'infix' is set to either 'model.model_type'
%                      if existent or the empty string.
%
% generated variables:
%   - 'detailed_data': structure containing the reduced basis space and
%                      information on its construction process.
%

disp('constructing reduced basis')
tmp = load(fullfile(rbmatlabresult,CRBfname));
detailed_data = tmp.detailed_data;

if structcmp(model, tmp.model) ~= 1
  warning('model has changed since crb basis generation.');
end


detailed_data.W = model_data.W;
if ~isfield(params,'Mstrich')
  params.Mstrich = 0;
end
if params.Mstrich == 0
  model.enable_error_estimator = 0;
end


if ~isfield(params,'M')
  model.M       = cellfun(@(x)(size(x,2) - model.Mstrich), detailed_data.BM, 'UniformOutput', true)
                   - params.Mstrich;
else
  model.M       = params.M;
end
model.Mstrich   = params.Mstrich;

tic;
detailed_data   = rb_basis_generation(model, detailed_data);
t = toc;

detailed_data.RB_info.elapsed_time = t;
if isempty(detailedfname)
  if isfield(model, 'model_type')
    infix = model.model_type;
  else
    infix = '';
  end
  detailedfname = [model.name, infix, '_detailed.mat'];
end
save(fullfile(rbmatlabresult,detailedfname),...
     'detailed_data','model','plot_params');

plot(detailed_data.RB_info.max_err_sequence);
set(gca,'Yscale','log');
title('RB-generation error convergence');

