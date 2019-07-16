% script comparing time for a reduced and a detailed simulation and
% starting the demonstration GUI
%
% required variables that need to be set:
%    - 'detailedfname'

tmp=load(fullfile(rbmatlabresult,detailedfname));
detailed_data=tmp.detailed_data;
disp('reduced simulation:')
reduced_data = gen_reduced_data(model, detailed_data);
model.N     = size(detailed_data.RB,2);
model.Mstrich = 0;
model.M     = cellfun(@(x)(size(x,2) - model.Mstrich), detailed_data.BM, 'UniformOutput', true);
model.enable_error_estimator = 0;
reduced_data = extract_reduced_data_subset(model, reduced_data);
tic;
rb_sim_data = rb_simulation(model, reduced_data);
t = toc;
disp(['time for online phase: t = ',num2str(t)]);

disp('full simulation:')
tic;
sim_data = detailed_simulation(model, detailed_data);
t = toc;
disp(['time for detailed simulation: t = ',num2str(t)]);

demo_rb_gui(model, detailed_data, [], plot_params);

