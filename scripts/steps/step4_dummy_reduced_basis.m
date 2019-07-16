% script constructing a reduced basis from a single trajectory and
% performing a detailed simulation where the data is projected on the
% dummy reduced basis space after each time step.
%
% required variables that need to be set:
%   - 'CRBfname'
%   - 'plot_params'
%
% generated variables:
%   - 'ei_rb_proj_data': output of detailed_ei_rb_proj_simulation() with
%                        dummy reduced basis space.

tmp=load(fullfile(rbmatlabresult,CRBfname));
detailed_data=tmp.detailed_data;
disp('detailed interpolated simulation for basis construction:')
if ~isfield(model,'Mstrich')
  model.Mstrich = 0;
  model.enable_error_estimator = 0;
end;

model.M = cellfun(@(x)(size(x,2) - model.Mstrich), detailed_data.BM, 'UniformOutput', true);

tic
ei_sim_data =  model.detailed_ei_simulation(model,detailed_data);
dtime = toc

% in case of newton scheme add temporary results of newton steps to dummy
% reduced basis
if model.newton_solver
  U = ei_sim_data.Unewton;
else
  U = ei_sim_data.U;
end

UON = orth(U);
UON = UON / (full(sqrt(detailed_data.W(1))));

%UON = model.orthonormalize(model, detailed_data, ei_sim_data.U,1e-7);
%
detailed_data.RB = UON;

disp('reduced simulation:')
reduced_data = gen_reduced_data(model, detailed_data);
model.N = size(detailed_data.RB,2);
smodel = model;
smodel.enable_error_estimator = 0;

tic
rb_sim_data = rb_simulation(smodel, reduced_data);
rbtime = toc
tic
rb_sim_data = rb_reconstruction(smodel, detailed_data, rb_sim_data);
rectime = toc

if ~isfield(rb_sim_data, 'Delta')
  smodel.enable_error_estimator = 0;
else
  Delta        = rb_sim_data.Delta(end);

  disp(['Error estimator returned: ', num2str(Delta)]);
end

disp('detailed interpolated and rb-projected simulation:')
model.M = cellfun(@(x)(size(x,2) - model.Mstrich), detailed_data.BM, 'UniformOutput', true);
model.N = size(detailed_data.RB,2);
%ei_rb_proj_data =  detailed_ei_rb_proj_simulation(smodel, detailed_data);

plot_params.title = 'reduced simulation result';
plot_sim_data(model, detailed_data, rb_sim_data, plot_params);
plot_params.title = 'detailed ei_interpol simulation result';
plot_sim_data(model, detailed_data, ei_sim_data, plot_params);
%plot_params.title = 'detailed ei_interpol and rb_projected simulation result';
%plot_sim_data(model, detailed_data, ei_rb_proj_data, plot_params);

epp = plot_params;
epp.title = 'error';
diff_data         = rb_sim_data;
%diff_data.U       = abs(diff_data.U - ei_rb_proj_data.U);
diff_data.U       = abs(diff_data.U - ei_sim_data.U);
epp.clim = [0,max(max(diff_data.U))];
plot_sim_data(model, detailed_data, diff_data, epp);
disp(['maximum l-infty error:',num2str(max(max(diff_data.U)))]);

if model.debug && smodel.enable_error_estimator
  ei_res       = rb_sim_data.eilog(end);
  max_ei_error = max(cellfun(@(X) X.max_err_sequence(model.M), ...
                           detailed_data.ei_info));
  ei_rel_diff  = abs(ei_res - max_ei_error) / abs(ei_res);

  if ei_rel_diff < 0.1 || ei_rel_diff > 10
    disp(['Oh no! The EI residuum has a different order of magnitude ', ...
          'than the empirical interpolation error.']);
  else
    disp(['EI residuum has same order of magnitude as empirical ', ...
          'interpolation error.']);
  end;
  disp(['EI residuum returned: ', num2str(ei_res)]);
  disp(['maximum EI error:     ', num2str(max_ei_error)]);

  proj_res      = rb_sim_data.reslog(end);
%  proj_res_exa  = ei_rb_proj_data.reslog(end);
%  proj_res_diff = abs(proj_res_exa - proj_res) / abs(proj_res);

%  if proj_res_diff < 0.1 || proj_res_diff > 10
%    disp(['Oh no! The exact projection residuum has a different order of ', ...
%          'magnitude than the estimated.']);
%  else
%    disp(['Exact projection residuum has same order of magnitude as ', ...
%          'estimated projection residuum.']);
%  end;
  disp(['estimated projection residuum: ', num2str(proj_res)]);
  disp(['exact projection residuum:     ', num2str(proj_res_exa)]);
end

disp(['runtime detailed       : ', num2str(dtime)]);
disp(['runtime rb_sim         : ', num2str(rbtime),  '  factor: ', num2str(dtime/rbtime)]);
disp(['runtime reconstruction : ', num2str(rectime), '  factor: ', num2str(dtime/(rbtime+rectime))]);

