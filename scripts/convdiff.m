% small script demonstrating the convdiff example from the M2AN Paper, that is
% also implemented in Dune. Later it will be possible to use the Dune
% implementation through a mex interface.

% Martin Drohmann 06.05.2009
% based on burgers_fv.m by
% Bernard Haasdonk 14.8.2007

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Select here, what is to be performed                        %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
step = 0 % initialize the model_data structure. This step needs to be executed
          % once per Matlab session only
%step = 1 % single detailed simulation with given data and plot. Run
          % this with varying parameters mu until sure that scheme
          % is stable. Modify dt or the data-functions accordingly,
          % until a nice parameter-domain with uniformly stable
          % detailed scheme is obtained.
%step = 2 % generate dummy reduced basis from single trajectory and check, if
          % ei_interpolation with projection on this space maintains
          % result. A simple reduced simulation can also be
          % performed. All results should be visually identical
%step = 3 % generate reduced basis
%step = 6 % time measurement of reduced simulation and
          % use reduced basis in rb_demo_gui
%step = 7 % generate error-landscape over varying N and M
          % can take several hours!!!
%step = 8 % do runtime comparisons between detailed and reduced simulations

%steps = {2,5,7,8}
steps = {0,1};

for si=1:length(steps)

  step = steps{si};

  % output-filenames in rbmatlabresult
  detailedfname = 'convdiff_detailed_interpol.mat';

  %% parameters for visualization
  plot_params.show_colorbar = 1;
  plot_params.colorbar_mode = 'EastOutside';
  plot_params.plot          = @fv_plot;

  switch step
   case 0 % initialize model data
    model            = convdiff_model;

    model_data = gen_model_data(model);
   case 1 % single detailed simulation and plot
    disp('performing single detailed simulation')
    tic
    mu_test   = cellfun(@mean, model.mu_ranges);
    sim_data  = detailed_simulation(model, model_data);
    toc
    plot_sim_data(model, model_data, sim_data, plot_params);

   case 2 % construct dummy reduced basis by single trajectory and simulate
    disp('detailed interpolated simulation for basis construction:')
    mu_test   = cellfun(@mean, model.mu_ranges);
    %  mu_test2  = cellfun(@min, params.mu_ranges);
    model = model.set_mu(model, mu_test);
    sim_data  = detailed_simulation(model, model_data);
    UON       = model.orthonormalize(model, model_data, sim_data.U);
    detailed_data.grid = model_data.grid;
    detailed_data.W    = model_data.W;
    detailed_data.RB = UON;
    disp('reduced simulation:')
    reduced_data = gen_reduced_data(model, detailed_data);
    model.N = size(detailed_data.RB,2);
    reduced_data = extract_reduced_data_subset(model, reduced_data);
    model = model.set_mu(model, mu_test);
    rb_sim_data = rb_simulation(model, reduced_data);
    rb_sim_data = rb_reconstruction(model, detailed_data, rb_sim_data);

    plot_params.title = 'reduced simulation result';
    plot_sim_data(model, model_data, rb_sim_data, plot_params);
    plot_params.title = 'detailed simulation result';
    plot_sim_data(model, model_data, sim_data, plot_params);

   case 3 % reduced basis
    mu_test   = cellfun(@mean, model.mu_ranges);
    disp('constructing reduced basis')
    detailed_data = gen_detailed_data(model, model_data);
    tic;
    %  detailed_data = rb_basis_generation(detailed_data, ...
    %                  params);
    t = toc;
    detailed_data.RB_info.elapsed_time = t;
    save(fullfile(rbmatlabresult,detailedfname),...
         'detailed_data','params');
    plot(detailed_data.RB_info.max_err_sequence);
    set(gca,'Yscale','log');
    title('RB-generation error convergence');
   case 4
    mu_test   = cellfun(@mean, model.mu_ranges);
    load(fullfile(rbmatlabresult,detailedfname));
    disp('reduced simulation:')
    reduced_data = gen_reduced_data(model, detailed_data);
    model.N = size(detailed_data.RB,2);

    reduced_data = extract_reduced_data_subset(model, reduced_data);
    tic;

    model = model.set_mu(model, mu_test);
    rb_sim_data = rb_simulation(model,reduced_data);
    t = toc;
    disp(['time for online phase: t = ',num2str(t)]);

    disp('full simulation:')
    tic;
    sim_data = detailed_simulation(model,model_data);
    t = toc;
    disp(['time for detailed simulation: t = ',num2str(t)]);

    demo_rb_gui(model,detailed_data,[],plot_params);

   case 7 % training-error landscape
    disp('warning: takes a few hours!');
    load(fullfile(rbmatlabresult,detailedfname));

    model.N = params.RB_stop_Nmax;

    offline_data = gen_reduced_data(model,detailed_data);
    range_params.plot_fields = { 'N', 'M' };
    range_params.max         = [size(detailed_data.RB,2), ...
                                size(detailed_data.BM{1},2) ];
    range_params.sample_size = [ 7, 12 ];
    range_params.mu_set_size = 100;

    testdir = ['convdiff_test_', range_params.mu_set_size ];

    params.run_name    = [testdir, '_step7'];
    params.tictoctable = true;

    output = stochastic_error_estimation(model, detailed_data, offline_data,...
                                         testdir, range_params, params);

  otherwise
    error('step-number is unknown!');

  end; % switch


end; % for si
%| \docupdate 
