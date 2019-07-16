function detailed_data = basisgen_refined_tmp(detailed_data, params)
%function detailed_data = basisgen_refined_tmp(detailed_data,params)
%
% script demonstrating the generation of a reduced basis with
% an refined parameter grid approach, either uniformly refined or 
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
% basis-extension again started. The new training set is only the
% set of new nodes.
%
% produced fields of detailed_data:
%           RB : collection of orthonormal reduced basis DOF
%                vectors
%           RB_info: datastructure with basis-generation
%                    information, see below
%           grid : numerical grid, is constructed if not existent.
%
% generated fields in RB_info in addition to those by RB_greedy_extension:
%           M_train_list : a cell array containing the training sets
%                          provided to RB_greedy_extension.
%           MMesh_list : a cell array containing the sequence of
%                        increasingly refined meshes   
%           mesh_index_sequence : a list indicating, on which refinement-level
%                                 each extension-mu-vector in 
%                                 mu_sequence was added.
%           RB_stopped_on_max_refinement_level : indicates, whether
%                            the maximum level of refinements was reached 
%
% required fields of params in addition to those of RB_greedy_extension:
%         RB_refinement_mode: 'uniform' or 'adaptive'
%         RB_max_refinement_level : level after which
%                  no further refinement is performed 
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

% generate grid
if ~isfield(detailed_data,'grid') | isempty(detailed_data.grid)
  detailed_data.grid = construct_grid(params); 
end;

if ~isfield(params,'RB_max_refinement_level')
  params.RB_max_refinement_level = inf; 
end;

%%%%%%%%%%%%%%% specific settings for adaptive algorithm

par.range = params.mu_ranges;
par.numintervals = [params.RB_numintervals];
MMesh0 = cubegrid(par);
Mfine = get(MMesh0,'vertex')';
Mcoarse = zeros(size(Mfine,1),0); % i.e. empty, but row-number as Mfine

% generate RB-start: all initial data constellations
detailed_data.RB = RB_init_data_basis(Mfine, ...
				      detailed_data.grid,...
				      params);

%prepare params 
params.Nmax = size(detailed_data.RB,2);

disp(['Starting RB extension loop']);

% start time measurement
tic;

% generate RB_info if required
if ~isfield(detailed_data,'RB_info')
  detailed_data.RB_info = [];
end;

% generate validation data if required
if ~isfield(detailed_data.RB_info,'M_val')
  rand('state',params.RB_val_rand_seed); % random but deterministic
  detailed_data.RB_info.M_val = rand_uniform(params.RB_M_val_size,...
					     params.mu_ranges);
end;

ref_loop = 0; % number of refinement loops
MMesh = MMesh0;
% monitor certain quantities during extension
MMesh_list = {};
M_train_list = {};
mesh_index_sequence = [];

% loop over refinement levels
continue_loop = 1;

while continue_loop

  % 1. Greedy search with current grid-vertices minus coarse grid-vertices
      
  % temporarily set M_train
  % M_fine = get(MMesh, 'vertex')'; 
  
  %  [M_train, i_fine_without_coarse] = vectorset_difference(Mfine,Mcoarse);
  M_train = Mfine;
  i_fine_without_coarse = 1:size(M_train,2);
  
  detailed_data.RB_info.M_train = M_train;
  
  disp(['Size of training set in next greedy: ',num2str(size(M_train,2))]);
%  disp('check set difference!')
%  keyboard;
  
  % extend basis
  [detailed_data,offline_data] = RB_greedy_extension(detailed_data,params);
  
  % extend mesh-information if required
  nRB_new = size(detailed_data.RB,2);
  lmi = length(mesh_index_sequence);
  if (nRB_new> lmi)
    mesh_index_sequence = [mesh_index_sequence, ...
		    (ref_loop+1) * ones(1,nRB_new-lmi)];
    MMesh_list = [MMesh_list, {MMesh}];
    M_train_list = [M_train_list, {M_train}];
  end;
  
  % determine eps max on Mfine cut Mcoarse  
  i_mask = zeros(1,size(Mfine,2));
  i_mask(i_fine_without_coarse) = 1;
  i_fine_cut_coarse = find(i_mask == 0);
  epsmax_fine_cut_coarse = - inf; 
  if ~isempty(i_fine_cut_coarse)
    Mfine_cut_coarse = Mfine(:,i_fine_cut_coarse);
    tmpparams = params;
    tmpparams.Nmax = size(detailed_data.RB,2);
    tmpparams.N = params.Nmax;
    errs_fine_cut_coarse = ...
	rb_test_indicator(detailed_data, offline_data,...
			  Mfine_cut_coarse,...
			  [],tmpparams);
    [max_errs, mu_inds] = max(errs_fine_cut_coarse);
    epsmax_fine_cut_coarse = max_errs(1);
  end;
  
  % construct eps_max on Mfine by both eps_max over (Mfine cut Mcoarse) and 
  % (Mfine \ Mcoarse))
  
  if isempty(detailed_data.RB_info.max_err_sequence)    
    % if no extension has been performed: get error from M_first_errs
    epsmax_fine_without_coarse = ...
	max(detailed_data.RB_info.M_first_errs);
  else % an extension has been performed    
    epsmax_fine_without_coarse = ...
	detailed_data.RB_info.max_err_sequence(end);
  end;
  
  epsmax_fine = max(epsmax_fine_cut_coarse, ...
		     epsmax_fine_without_coarse);

  continue_loop = 1;
  
  %  check if error on current grid is sufficient
  stopped_on_epsilon = 0;
  if epsmax_fine < params.RB_stop_epsilon
    stopped_on_epsilon = 1;
    continue_loop = 0;
  end;
  
  stopped_on_max_refinement_level = 0;
  levels = get(MMesh,'level');
  maxlev = max(levels);
  if maxlev >= params.RB_max_refinement_level
    stopped_on_max_refinement_level = 1;
  end;
  
  skip_refinement_and_break = stopped_on_epsilon | ...
      detailed_data.RB_info.stopped_on_timeout | ...
      detailed_data.RB_info.stopped_on_Nmax | ...
      stopped_on_max_refinement_level;
  
  if skip_refinement_and_break
    continue_loop = 0;
  else   
    % 2. Grid refinement
        
    disp('detected grid refinement necessary!!');    

    nleafelements = get(MMesh,'nleafelements');
    
    switch params.RB_refinement_mode
     case 'uniform'
      
      MMesh = refine(MMesh,1:nleafelements);   
      
     case 'adaptive'
      
      % compute error estimators       
      M_last_errs = -1 * ones(1,size(Mfine,2));
      % last call of RB_greedy_extension gave fine_without_coarse
      M_last_errs(i_fine_without_coarse) = detailed_data.RB_info.M_last_errs;
      % call of rb_test_indicator gave fine_cut_coarse
      if ~isempty(i_fine_cut_coarse)
	M_last_errs(i_fine_cut_coarse) = errs_fine_cut_coarse;
      end;
      
      if ~isempty(find(M_last_errs<0))
	error('errors for refinement-estimators not all specified!!');
      end;
      
      %for all leaf e in grid M_t compute
      %     eta(e) := 
      %         max_{mu in (V(e) and cog(e)   Delta (Phi_t, mu) 
      %determine Theta most violating elements

      % eta: pol0-function on mu-grid-elements, only leaf elements
      [eta, eta_info] = rb_mu_element_error_estimators ( ...
	  detailed_data, offline_data, MMesh, ...
	  M_last_errs, ...
	  params);
      
      disp(['saving before grid-refinement ',num2str(ref_loop)]);
      fn = ['tmp_adapt',num2str(ref_loop),'.mat'];
      save(fullfile(getenv('PDEMATLABTEMP'),fn));
      
      [dummy , ind] = sort(-eta);
      nmax = ceil(nleafelements * params.RB_refinement_theta);
      li = get(MMesh,'leafelements');
      MMesh = refine(MMesh,ind(1:nmax));   
      
     otherwise
      error('refinement mode unknown!');
    end;      

    % compute new Mfine 
    
    Mcoarse = Mfine;
    Mfine = get(MMesh, 'vertex')'; 
    
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
detailed_data.RB_info.M_train_list  = M_train_list;
detailed_data.RB_info.MMesh_list = MMesh_list;

detailed_data.RB_info.RB_stopped_on_max_refinement_level = ... 
    stopped_on_max_refinement_level;

detailed_data.RB_info.RB_stopped_on_epsilon = ... 
    stopped_on_epsilon;

disp(['Generated RB basis on refined grid with ',...
      num2str(size(detailed_data.RB,2)), ...
      ' basis vectors.']);
disp(['Number of MMesh-refinement loops = ',num2str(ref_loop)]);
t = toc;
disp(['Runtime = ',num2str(t)]);


 




% TO BE ADJUSTED TO NEW SYNTAX
%| \docupdate 
