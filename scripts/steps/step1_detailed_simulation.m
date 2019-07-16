% script performing a single detailed simulation and plotting it.
%
% required variables that need to be set:
%   - 'model'
%   - 'model_data'
%   - 'plot_params'
%
% optional variables that can be set:
%   - 'mu_default': a mu vector for which the detailed simulation shall
%                   be computed.
%
% generated variables:
%   - 'sim_data'

disp('performing single detailed simulation')

%load(fullfile(rbmatlabtemp, params.step0file));

tic;
%if isempty(mu_default)
%  mu_default = cellfun(@mean, model.mu_ranges);
%end
%model = model.set_mu(model, mu_default);
sim_data = detailed_simulation(model, model_data);
toc
plot_sim_data(model, model_data, sim_data, plot_params);

%save(fullfile(rbmatlabtemp, params.step1file), 'sim_data');

