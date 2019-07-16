function detailed_data = rb_basis_generation(model,detailed_data)
%function detailed_data = rb_basis_generation(model,detailed_data)
% reduced basis construction with different methods
%
% reduced basis construction with different methods. Performs
% computing and insertion of a field RB and RB_info into the
% detailed_data structure. detailed_data is assumed to contain the
% grid if required
%
% Required fields of model:
%  mu_names : cell array of parameter-names
%  mu_ranges : cell array of parameter-intervals
%  RB_generation_mode :
%  - \c 'file_load' : simple loading of reduced basis from a file
%  - \c 'lagrangian' : simple snapshots without orthonormalization
%  - \c 'PCA_trajectory' : simple detailed simulation and PCA of the
%  - \c 'PCA_trajectories' : simple detailed simulation and PCA of the
%                     trajectories with parameter vectors specified in
%                     model.RB_mu_list
%  - \c 'greedy_uniform_fixed' : greedy algorithm based on a uniform
%                           cartesian parameter grid
%  - \c 'greedy_log_uniform_fixed' : greedy algorithm based on a
%                           logarithmically uniform
%                           cartesian parameter grid
%  - \c 'greedy_random_fixed' :  greedy algorithm based on a uniformly
%                           distributed random parameter set
%  - \c 'greedy_refined' :  greedy algorithm based on a uniformly
%                       random parameter set with grid refinement
%  - \c 'model_RB_basisgen' : call individual model.RB_basisgen 
%  - \c 'none' : skipping of basis generation
%
% Optional fields of model:
%    RB_basis_filename : name of file, that contains precomputed
%                        reduced basis RB in 'file_load' mode
%    RB_train_rand_seed : set random seed for M_train generation in random mode
%    RB_train_size : number of training parameter vectors to ...
%                    generate in random-mode
%    RB_numintervals: vector indicating the number of intervals of
%                     the (coarse) grid in parameter space
%
% generated fields of detailed_data:
%  RB    : reduced basis columns 1,...,Nmax
%  RB_info : depending on generation method some detailed information

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


% Bernard Haasdonk 15.5.2007

% import necessary model methods

% initialize empty, such that get_rb_size works subsequently
detailed_data = model.set_rb_in_detailed_data(detailed_data,[]);

switch model.RB_generation_mode
 case 'file_load'

  tmp = load(model.RB_basis_filename);
  detailed_data.RB_info = [];
  detailed_data = model.set_rb_in_detailed_data(detailed_data,tmp.RB);
%  detailed_data.RB = tmp.RB;

 case 'lagrangian'
  
  Utot = [];
  for m = 1:length(model.RB_mu_list)
    model = set_mu(model,model.RB_mu_list{m});
    sim_data = detailed_simulation(model,detailed_data);
    if isempty(Utot)
      Utot = model.get_dofs_from_sim_data(sim_data);
    else
      U = model.get_dofs_from_sim_data(sim_data);
      Utot = [Utot U];
    end;
  end;
  detailed_data = model.set_rb_in_detailed_data(detailed_data,Utot);
  
 case 'from_detailed_data'

  detailed_data.RB_info = [];
  detailed_data = model.set_rb_in_detailed_data(...
      detailed_data,detailed_data.RB);
  %  detailed_data.RB = tmp.RB;

 case 'PCA_trajectory'

  detailed_data.RB_info = [];
  sim_data = detailed_simulation(model,detailed_data);
  if ~isfield(model,'RB_stop_Nmax')
    RB = model.PCA_fixspace(model.get_dofs_from_sim_data(sim_data), ...
                            []);
  else
    RB = model.PCA_fixspace( ...
        model.get_dofs_from_sim_data(sim_data), ...
        [],...
        model.get_inner_product_matrix(detailed_data), ...
        model.RB_stop_Nmax);
  end;
  detailed_data = model.set_rb_in_detailed_data(detailed_data,RB);

 case 'PCA_trajectories'

  Utot = [];
  for m = 1:length(model.RB_mu_list)
    model = set_mu(model,model.RB_mu_list{m});
    sim_data = detailed_simulation(model,detailed_data);
    if isempty(Utot)
      Utot = model.get_dofs_from_sim_data(sim_data);
    else
      U = model.get_dofs_from_sim_data(sim_data);
      Utot = [Utot U];
    end;
  end;
  detailed_data.RB_info = [];
  W = model.get_inner_product_matrix(detailed_data);
  %%%  gram_schmidt too inaccurate in the following!!!
  %  RB = orthonormalize_qr(Utot,W);
  %  disp('perhaps also try PCA_fixspace in this case!!');

  if ~isfield(model,'RB_stop_Nmax')
    RB = model.PCA_fixspace(Utot, ...
                            [],W,[],'qr');
  else
    RB = model.PCA_fixspace( ...
        Utot, ...
        [],...
        W, ...
        model.RB_stop_Nmax,'qr');
  end;

  detailed_data=model.set_rb_in_detailed_data(detailed_data,RB);

 case 'greedy_random_fixed'

  rand('state',model.RB_train_rand_seed);
  %model.RB_extension_stop = 'epsilon';
  detailed_data.RB_info = [];
  detailed_data.RB_info.M_train = rand_uniform(model.RB_train_size,...
                                                  model.mu_ranges);
  detailed_data = basisgen_fixed(model,detailed_data);

 case 'greedy_uniform_fixed'

  par.numintervals = model.RB_numintervals;
  par.range = model.mu_ranges;
  %  enable restart of basis generation!!!
  %  detailed_data.RB_info = [];
  
  MMesh0 = cubegrid(par);
  detailed_data.RB_info.M_train = get(MMesh0,'vertex')';
  %detailed_data.RB_info.M_train = [0;1;0];
  detailed_data = basisgen_fixed(model,detailed_data);

 case 'greedy_log_uniform_fixed'

  par.numintervals = model.RB_numintervals;
  par.range = model.mu_ranges;
  for i = 1:length(par.range);
    par.range{i} = log(par.range{i});
  end;
  %  enable restart of basis generation!!!
  %  detailed_data.RB_info = [];
  MMesh0 = cubegrid(par);
  Mtrain_log = get(MMesh0,'vertex')';
  detailed_data.RB_info.M_train = exp(Mtrain_log);
%  keyboard; % check
  detailed_data = basisgen_fixed(model,detailed_data);

 case 'greedy_refined'

  detailed_data.RB_info = [];
  par.numintervals = model.RB_numintervals;
  par.range = model.mu_ranges;
  %  enable restart of basis generation!!!
  %  detailed_data.RB_info = [];
  MMesh0 = cubegrid(par);
  %  detailed_data.RB_info.M_train = get(MMesh0,'vertex')';
  detailed_data.RB_info.MMesh0 = MMesh0;
  detailed_data = basisgen_refined(model, detailed_data);
%  detailed_data.RB_info.M_val = rand_uniform(model.RB_M_val_size,...
%                                            model.mu_ranges);

 case 'model_RB_basisgen'
  detailed_data.RB = model.RB_basisgen(model,detailed_data);
  
 case 'none'
  disp(['warning: no reduced basis is produced, detailed_data will' ...
        ' not be suitable for simulation now!!']);

 otherwise
  error('RB_generation_mode unknown!!!');
end;

if isfield(model,'debug')
if model.debug && isfield(detailed_data, 'RB') && ~isempty(detailed_data.RB)
  K = detailed_data.RB'* ...
      model.get_inner_product_matrix(detailed_data) * detailed_data.RB;
  e = max(max(abs(K-eye(size(K)))));
  if e>1e-4
    error('reduced basis not orthonormal!!')
  end;
    end;
end;
