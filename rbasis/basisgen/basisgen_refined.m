function detailed_data = basisgen_refined(model,detailed_data)
%function detailed_data = basisgen_refined(model,detailed_data)
%
% script demonstrating the generation of a reduced basis with
% a refined parameter grid approach, either uniformly refined or 
% adaptively.
%
% RB is started with initial data by varying all parameters as
% specified in the coarse parameter grid M_0.
% For current parameter mesh M_t, the basis is extended by
% greedy searching for the node with the maximum posterior error
% estimator, for this mu a single vector is chosen for basisextension. 
% A validation error is determined for checking, whether overfitting is
% currently happening. If so, the grid is refined either uniformly
% or by an a-posteriori error estimation approach, and the
% basis-extension again started.
%
% produced fields of detailed_data:
%           RB,V,W, etc. : collection of orthonormal reduced basis DOF
%                vectors
%           RB_info: datastructure with basis-generation
%                    information, see below
%           grid : numerical grid, is constructed if not existent.
%
% generated fields in RB_info:
%           M_train : is deleted and only set temporarilty in this function
%           MMesh_list : a cell array containing the sequence of
%                        increasingly refined meshes   
%           mesh_index_sequence : a list indicating, on which refinement-level
%                                 each extension-mu-vector in 
%                                 mu_sequence was added.
%           RB_stopped_on_max_refinement_level : indicates, whether
%                            the maximum level of refinements was reached 
%           RB_max_val_train_ratio_assured : indicates whether the maximum
%                     validation training ratio has been fulfilled at any time.
%
% required fields of model:
%         RB_refinement_mode: 'uniform' or 'adaptive'
%         RB_max_refinement_level : level after which
%                  no further refinement is performed (Default: 'inf')
%         RB_force_stop_max_refinement_level : flag indicating whether the
%                  basis generation algorithm needs to be stopped when the
%                  maximum refinment level is reached. If set to 'false',
%                  'RB_stop_max_val_train_ratio' is set to '1e6' after last
%                  possible refinement. (Default: false)
%         mu_ranges : cell array of intervals for the mu-components
%         RB_numintervals : vector indicating the number of
%                       intervals for the coarse initial mesh
%         RB_val_rand_seed: random seed for validation set
%                      creation, if not already given in RB_info.M_val
%         RB_M_val_size: number of validation points for extension
%         RB_refinement_theta: theta-value in [0,1] in case of
%                              adaptive refinement

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


% Bernard Haasdonk 29.5.2007


% We assume, that detailed_data is reasonably filled, e.g. contains
% the grid.

% generate grid
%if ~isfield(detailed_data,'grid') | isempty(detailed_data.grid)
%  detailed_data.grid = construct_grid(model); 
%end;

if isfield(model,'RB_stop_max_refinement_level')
    error('field RB_stop_max_refinement_level has been replaced by RB_max_refinement_level. If you want to force RB generation to stop use RB_force_stop_max_refinement_level')
end

if ~isfield(model,'RB_max_refinement_level')
  model.RB_max_refinement_level = inf; 
end;

if ~isfield(model,'RB_force_stop_max_refinement_level')
  model.RB_force_stop_max_refinement_level = false;
end

%%%%%%%%%%%%%%% specific settings for adaptive algorithm

%par.range = model.mu_ranges;
%par.numintervals = [model.RB_numintervals];
%MMesh0 = cubegrid(par);
MMesh0 = detailed_data.RB_info.MMesh0;
detailed_data.RB_info.M_train = get(MMesh0,'vertex')';

% generate RB-start: all initial data constellations

if (model.get_rb_size(model,detailed_data)==0)  
  % generate RBstart: all initial data constellations
  RB = model.rb_init_data_basis(model, detailed_data);
  detailed_data = model.set_rb_in_detailed_data(detailed_data,RB);
  %  if (get_rb_size(detailed_data)==0) % check if reduced basis is empty
  %  detailed_data = model.RB_init_data_basis(model, detailed_data, ...
  %                                           detailed_data.RB_info.M_train);
end;

%detailed_data.RB = RB_init_data_basis(M_train, ...
%				      detailed_data.grid,...
%				      model);


%prepare model 

%model.Nmax = size(detailed_data.RB,2);
model.Nmax = model.get_rb_size(model,detailed_data);

disp('Starting RB extension loop');

% start time measurement
tic;
start_time_basisgen_refined = tic;

% generate RB_info if required
if ~isfield(detailed_data,'RB_info')
  detailed_data.RB_info = [];
end;

% generate validation data if required
if ~isfield(detailed_data.RB_info,'M_val')
  rand('state',model.RB_val_rand_seed); % random but deterministic
  detailed_data.RB_info.M_val = rand_uniform(model.RB_M_val_size,...
					     model.mu_ranges);
end;

ref_loop = 0; % number of refinement loops
MMesh = MMesh0;
% monitor certain quantities during extension
MMesh_list = {};
mesh_index_sequence = [];

% loop over refinement levels
continue_loop = 1;

max_val_train_ratio = model.RB_stop_max_val_train_ratio;

while continue_loop

  % 1. Greedy search with current grid-vertices
  
  % temporarily set M_train
  % M_fine = get(MMesh, 'vertex')'; 
  
  detailed_data.RB_info.M_train = get(MMesh, 'vertex')'; 
  
  % extend basis
  [detailed_data,offline_data] = RB_greedy_extension(model,detailed_data);

  % extend mesh-information if required
  nRB_new = model.get_rb_size(model,detailed_data);
  lmi = length(mesh_index_sequence);
  if (nRB_new > lmi)
    mesh_index_sequence = [mesh_index_sequence, ...
		    (ref_loop+1) * ones(1,nRB_new-lmi)];
    MMesh_list = [MMesh_list, {MMesh}];
  end;

  stopped_on_max_refinement_level = 0;
  levels = get(MMesh,'level');
  maxlev = max(levels);

  if maxlev == model.RB_max_refinement_level - 1 ...
      && ~model.RB_force_stop_max_refinement_level
    model.RB_stop_max_val_train_ratio = 1e6;
  end

  % force stop of basis generation after maximum refinement level is reached
  if maxlev >= model.RB_max_refinement_level ...
      && model.RB_force_stop_max_refinement_level
    stopped_on_max_refinement_level = 1;
  end

  % if refinement is still possible and required val-train ratio is not assured
  % refine the parameter space regardless of stopped_on_epsilon flag.
  if maxlev < model.RB_max_refinement_level ...
      && detailed_data.RB_info.r_value_sequence(end) > max_val_train_ratio
    detailed_data.RB_info.stopped_on_epsilon = 0;
  end;


  skip_refinement_and_break = detailed_data.RB_info.stopped_on_epsilon | ...
      detailed_data.RB_info.stopped_on_empty_extension | ...
      detailed_data.RB_info.stopped_on_timeout | ...
      detailed_data.RB_info.stopped_on_Nmax | ...
      stopped_on_max_refinement_level;
  
  if skip_refinement_and_break
    continue_loop = 0;
  else   
    % 2. Grid refinement
        
    disp('detected grid refinement necessary!!');    

    nleafelements = get(MMesh,'nleafelements');
    
    switch model.RB_refinement_mode
     case 'uniform'

      MMesh = refine(MMesh,1:nleafelements);   
      
     case 'adaptive'
      
      % compute error estimators       

      %for all leaf e in grid M_t compute
      %     eta(e) := 
      %         max_{mu in (V(e) and cog(e)   Delta (Phi_t, mu) 
      %determine Theta most violating elements

      % eta: pol0-function on mu-grid-elements, only leaf elements
      [eta, eta_info] = rb_mu_element_indicators ( ...
	  detailed_data, offline_data, MMesh, ...
	  detailed_data.RB_info.M_last_errs, model);
      
      %disp(['saving before grid-refinement ',num2str(ref_loop)]);
      %fn = ['tmp_adapt',num2str(ref_loop),'.mat'];
      %save(fullfile(getenv('PDEMATLABTEMP'),fn));
      
      [dummy , ind] = sort(-eta);
      nmax = ceil(nleafelements * model.RB_refinement_theta);
      %li = get(MMesh,'leafelements');
      MMesh = refine(MMesh,ind(1:nmax));   

     otherwise
      error('refinement mode unknown!');
    end;      
    
    ref_loop = ref_loop + 1;
    
    disp('finished grid refinement!!');
    
  end % end skip_refinement_and_break select
  
end % end loop over refinement levels



% prepare output:
detailed_data.RB_info = rmfield(detailed_data.RB_info,'M_train');

% the following field does also not make much sense in case of
% multiple greedy runs, so remove
detailed_data.RB_info = rmfield(detailed_data.RB_info,'M_first_errs');

detailed_data.RB_info.mesh_index_sequence  = mesh_index_sequence;
detailed_data.RB_info.MMesh_list = MMesh_list;

detailed_data.RB_info.RB_stopped_on_max_refinement_level = ... 
    stopped_on_max_refinement_level;
detailed_data.RB_info.RB_max_val_train_ratio_assured = ...
  (detailed_data.RB_info.r_value_sequence(end) < max_val_train_ratio);

disp(['Generated RB basis on refined grid with ',...
      num2str(model.get_rb_size(model,detailed_data)), ...
      ' basis vectors.']);
disp(['Number of MMesh-refinement loops = ',num2str(ref_loop)]);
t = toc(start_time_basisgen_refined);
disp(['Runtime = ',num2str(t)]);


%| \docupdate 
