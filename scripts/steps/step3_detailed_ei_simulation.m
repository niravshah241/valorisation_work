% script performing a detailed simulation with empirical interpolated
% operators comparing the result with a regular detailed simulation.
%
% required variables that need to be set:
%   - 'plot_params'
%   - 'CRBfname'
%   - 'mu_default'
%
% generated variables:
%   - 'sim_data':    output of regular detailed simulation
%   - 'ei_sim_data': output of detailed simulation with empirical
%                    interpolated operators.

disp(['comparison between detailed simulation with and without', ...
  ' interpolation']);
tmp = load(fullfile(rbmatlabresult,CRBfname));
detailed_data = tmp.detailed_data;

%model = model.set_mu(model, mu_default);
%sim_data = detailed_simulation(model, detailed_data);
plot_params.title = 'detailed simulation without interpolation';
plot_sim_data(model, detailed_data, sim_data, plot_params);
model.M = cellfun(@(x)(size(x,2) - model.Mstrich), detailed_data.BM, 'UniformOutput', true);

%model = model.set_mu(model, mu_default);
ei_sim_data =  detailed_ei_simulation(model, detailed_data);
plot_params.title = 'detailed simulation with empirical interpolation';
plot_sim_data(model, detailed_data, ei_sim_data, plot_params);
plot_params.title = 'error';
diff_data = sim_data;
diff_data.U = abs(diff_data.U - ei_sim_data.U);
diff_plot_params = plot_params;
diff_plot_params.clim = [0-eps, max(max(diff_data.U))];
plot_sim_data(model, detailed_data, diff_data, diff_plot_params);
disp(['maximum l-infty error:',num2str(max(max(diff_data.U)))]);

