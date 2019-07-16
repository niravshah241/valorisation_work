function detailed_data = nonlin_evol_gen_detailed_data(model, model_data)
%function detailed_data = nonlin_evol_gen_detailed_data(model, model_data)
% prepares detailed_data structure with high dimensional like reduced basis
% functions.
%
% method, which prepares 'detailed_data', which is meant as data, that may be
% dependent on `H`, and is not required during the online-simulation, but it
% may be used for reconstruction purposes or computation of online_quantities.
% i.e. reduced basis vectors, colateral reduced basis spaces, grid, etc. can be
% stored here.
%
% The grid, the collateral basis for the explicit nonlinearity and the reduced
% basis are generated in this order.
%
% \note The reduced basis generation is delegated to the function
% rb_basis_generation().
%
% - allowed dependency of data: `H`
% - allowed dependency of computation: `H`
% - unknown at this stage: `\mu, N_{\max}, M_{\max}, N, M`
%
% Required fields of model:
%  Nmax                : maximum number of reduced basis vectors to generate
%  Mmax                : maximum number of collateral reduced basis vectors to
%     generate
%  Mstrich             : number of collateral reduced basis function used for
%     a posteriori error estimation.
%  ei_numintervals     : (only used in 'param-time-space-grid' mode) the finite
%     parameter set for which the empirical interpolation is trained, is given
%     by through a uniform distribution of parameters in the bounds given by
%     'model.mu_ranges'. This field specifies the number of parameters in each
%     component direction.
%  ei_space_operators  : a cell array of function pointers to the localized
%     operators for which a collateral reduced basis space shall be created.
%  CRB_generation_mode : a string which is either
%     - 'file_load' : simple loading of collateral reduced basis from a file
%     - 'param-time-space-grid' : search for CRB vectors in a discrete subspace
%       of the space of all operator evaluations on solutions for parameters in
%       the 'model.mu_range' boundaries and all timesteps.
%  CRB_basis_filename  : (only used in 'file_load' mode) name of file, that
%      contains precomputed collateral reduced basis 'QM'.
%  separate_CRBs       : boolean flag, indicating whether, the operator
%      evaluations shall be stored in different directories for the later
%      computation of separate collateral reduced basis spaces for each
%      operator, or combined a single directory.
%
% Optional fields of model:
%  adaptive_time_split_mode : boolean flag indicating whether the time slices
%      shall be computed adaptively. (default == false)
%  time_split_Mmax          : when collateral reduced basis reaches this
%       size, the current time slice is split in 'adaptive_time_split_mode'
%       (default = 'model.ei_Mmax')
%  minimum_time_slice       : when time slice contains only this many number of
%       time steps, adaptation is stopped in 'adaptive_time_split_mode'.
%       (default = 1)
%  ei_stop_epsilon          : interpolation error limit stopping the basis
%       extension when reached. This needs to be strictly positive in
%       'adaptive_time_split_mode'. (default = 0 respectively 1e-6 in
%       'adaptive_time_split_mode')
%  ei_time_splits           : number of time slices the discrete time
%      interval shall (initially) be splitted into. (default = 1)
%
% generated fields of detailed_data:
%  grid  : grid to be used in the subsequent stages, copied from 'model_data'
%  QM    : a matrix of columnwise DOF-vectors `q_j` for interpolation
%          to be used as the colateral basis
%  TM    : a vector for the set of 'magic points' by containing the index
%          numbers of the corresponding DOF-nodes (In case of piecewise
%          constant or linear basis-functions, the maxima `x_i` always can be
%          found in a cell centroid (deg=0) or a node (deg=1).
%  BM    : the corresponding interpolation-matrix of dimension 'Mmax x Mmax',
%          i.e. interpolation can be done by solving the equation system `B
%          \sigma = \left(\zeta(x_1), \ldots, \zeta(x_M)\right)` then `Q
%          \sigma` is the empirical interpolation. Note, that this
%          reconstruction may not be done in the online phase!
%  ei_info: details on the empirical interpolation procedure as produced
%           by ei_detailed() plus the field 'extension_mus'
%  ei_info.extension_mus: a matrix whose column vectors hold the parameter
%           vectors corresponding to the solution trajectory serving as a basis
%           extension function.
%  RB    : matrix with 'Nmax' columns holding the reduced basis vectors
%
% Return values:
%  detailed_data : structure containing all the reduced basis data and
%           information about the generation process.

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

detailed_data = [];
detailed_data = structcpy(detailed_data, model_data);
params = [];

if ~isfield(model, 'adaptive_time_split_mode')
  model.adaptive_time_split_mode = false;
end

if ~isfield(model, 'ei_stop_epsilon')
  if model.adaptive_time_split_mode
    model.ei_stop_epsilon = 1e-6;
  else
    model.ei_stop_epsilon = 0;
  end
end

if ~isfield(model, 'ei_time_splits')
  model.ei_time_splits = 1;
end

if ~isfield(model, 'extend_crb')
  model.extend_crb = false;
end

if ~isfield(model, 'minimum_time_slice')
  model.minimum_time_slice = 1;
end

if ~isfield(model, 'time_split_Mmax')
  model.time_split_Mmax = model.Mmax;
end

if isa(model.ei_space_operators, 'function_handle')
    params.ei_space_operators = { model.ei_space_operators };

%if ~isfield(model, 'ei_space_operators')
%  params.ei_space_operators = { model.L_E_local_ptr };
%  detailed_data.explicit_crb_index = 1;
%  detailed_data.implicit_crb_index = -1;
%else
else
  params.ei_space_operators = model.ei_space_operators;
end
  if model.separate_CRBs
    detailed_data.implicit_crb_index = 2;
    detailed_data.explicit_crb_index = 1;
  else
    detailed_data.implicit_crb_index = 1;
    detailed_data.explicit_crb_index = 1;
  end
%end

if ~isfield(model, 'ei_time_splits')
  model.ei_time_splits = 1;
end

% empirical interpolation of nonlinearity must be performed before
% basis-construction, as basis construction needs availability of
% reduced simulation.

par.range = model.mu_ranges;

switch model.CRB_generation_mode

  % choose the following number of setting, that the time-snapshots
  % for the parameter ranges fit into memory:
  % 5x5x5 times 5 timesnapshots of size 8000 = 1250 x 8000 should be fine.
  % single nonlin simulation of 99 sec extrapolates to 3.5 hours for
  % 125 detailed simulations!

 case 'param-time-space-grid'

  disp('computing param-time-space-grid for empirical interpolation');
  par.numintervals = model.ei_numintervals;
  % if ei_operator_collect_tmp is used, set the following params!!
  %par.numintervals = [1,1,1]; 
  M = cubegrid(par);
  Mtrain = transpose(M.vertex);

  % initially set-up model.ei_time_splits time slices.
  detailed_data.time_split_map = [ (1:model.nt+1)', ...
    (ceil((1:model.nt+1)' / (model.nt+1)* (model.ei_time_splits-0.01))) ];

  % collect the operator evaluations
  LU_fnames = ei_operator_collect_files(model, detailed_data, Mtrain, ...
                                        params);

  ei_detailed_data = detailed_data;
  for lsi = 1:size(LU_fnames,1)
    params.ei_Mmax         = model.Mmax;
    if model.separate_CRBs
      disp(['computing CRB for local operator ', ...
            func2str(params.ei_space_operators{lsi})]);
      params.ei_Mmax = model.sMmax{lsi};
    else
      disp(cell2mat(['computing CRB for local operators: ', ...
            cellfun(@(X)([func2str(X),' ']), params.ei_space_operators, ...
                    'UniformOutput',false)]));
    end
    params.index     = lsi;
    if model.adaptive_time_split_mode
      max_time_splits = model.ei_time_splits;
      ti = 1;
      while ti <= max_time_splits;
        params.time_index      = ti;
        something_todo         = true;
        params.ei_Mmax         = model.time_split_Mmax;
        params.ei_stop_epsilon = model.ei_stop_epsilon;
        while something_todo
          ei_detailed_data  = ei_detailed(model, ei_detailed_data, ...
                                          LU_fnames(lsi,:), params);

          tei_info = ei_detailed_data.ei_info{lsi,ti};
          if tei_info.stopped_on_epsilon            % goto next time slice
            ti = ti + 1; something_todo = false;
            if model.verbose > 11
              slice_size = sum(ei_detailed_data.time_split_map(:,2) == ti);
              disp(['Reached epsilon for time slice ', num2str(ti), ...
                ' consisting of ' slice_size ' time steps.']);
            end
          elseif tei_info.stopped_on_Mmax % Have we reached Mmax?
            % yes, we reached the time_split_Mmax
            tmap = ei_detailed_data.time_split_map(:,2);
            slice_size = sum(tmap == ti);
            % do we have minimum slice size reached?
            if floor(slice_size/2) <= model.minimum_time_slice ...
                || tei_info.max_err_sequence(end) < model.ei_stop_epsilon_first % or first epsilon is reached
              if model.verbose > 11
                disp(['Reached minimum size for time slice ', num2str(ti), ...
                      ' consisting of ', num2str(slice_size), ' time steps.']);
              end
              % add further basis vectors, if requested...
              if model.time_split_Mmax ~= model.Mmax
                params.ei_Mmax = model.Mmax;
                ei_detailed_data  = ei_detailed(model, ei_detailed_data, ...
                                                LU_fnames(lsi,:), params);
              end
              % go to next time slice
              ti = ti + 1; something_todo = false;

            else
              % split the time slice!
              old_ti_indices = find(tmap == ti);
              minold = min(old_ti_indices);
              maxold = max(old_ti_indices);
              new_tmap = tmap;

              slmiddle = max(floor(median(tei_info.extension_filepos)), model.minimum_time_slice);
              slmiddle = min(slmiddle, slice_size-model.minimum_time_slice);
              %slmiddle = floor(slice_size/2);

              if model.verbose > 11
                disp(['Split time slice no. ', num2str(ti), ' at midpoint = ', ...
                      num2str(slmiddle) , '/', num2str(slice_size) ]);
              end


              new_tmap(minold:minold+slmiddle-1)   = ti;
              new_tmap(minold+slmiddle:maxold)     = ti+1;
              new_tmap(maxold+1:end) = tmap(maxold+1:end)+1;

              ei_detailed_data.time_split_map = [ (1:model.nt+1)', new_tmap ];

              disp(['new time split map: ', num2str(ei_detailed_data.time_split_map(:,2)')]);

              % find basis functions in first half and second half
              fh = tei_info.extension_filepos <= slmiddle;
              sh = ~fh;
              if model.extend_crb
                splitdd.QM = { ei_detailed_data.QM{lsi, ti}(:,fh), ...
                  ei_detailed_data.QM{lsi, ti}(:,sh) };
                splitdd.TM = { ei_detailed_data.TM{lsi, ti}(fh), ...
                  ei_detailed_data.TM{lsi, ti}(sh) };
                splitdd.BM = { ei_detailed_data.BM{lsi, ti}(fh,fh), ...
                  ei_detailed_data.BM{lsi, ti}(sh,sh) };
                splitdd.ei_info = cell(1,2);
                splitdd.ei_info{1,1}.max_err_sequence  = [tei_info.max_err_sequence(fh);Inf];
                splitdd.ei_info{1,2}.max_err_sequence  = [tei_info.max_err_sequence(sh);Inf];
                splitdd.ei_info{1,1}.extension_filenum = tei_info.extension_filenum(fh);
                splitdd.ei_info{1,2}.extension_filenum = tei_info.extension_filenum(sh);
                splitdd.ei_info{1,1}.extension_filepos = tei_info.extension_filepos(fh);
                splitdd.ei_info{1,2}.extension_filepos = tei_info.extension_filepos(sh)-slmiddle;
                params.skip_search = true;
              else
                splitdd.QM      = cell(1,2);
                splitdd.TM      = cell(1,2);
                splitdd.BM      = cell(1,2);
                splitdd.ei_info = cell(1,2);
                splitdd.ei_info{1,1}.max_err_sequence  = [];
                splitdd.ei_info{1,2}.max_err_sequence  = [];
                splitdd.ei_info{1,1}.extension_filenum = [];
                splitdd.ei_info{1,2}.extension_filenum = [];
                splitdd.ei_info{1,1}.extension_filepos = [];
                splitdd.ei_info{1,2}.extension_filepos = [];
              end

              newQMsize = size(ei_detailed_data.QM,2) + 1;
              ei_detailed_data.QM{lsi, newQMsize}      = [];
              ei_detailed_data.BM{lsi, newQMsize}      = [];
              ei_detailed_data.TM{lsi, newQMsize}      = [];
              ei_detailed_data.ei_info{lsi, newQMsize} = [];

              ei_detailed_data.QM(lsi, ti:end)         = ...
                {splitdd.QM{:}, ei_detailed_data.QM{lsi,ti+1:end-1}};
              ei_detailed_data.BM(lsi, ti:end)         = ...
                {splitdd.BM{:}, ei_detailed_data.BM{lsi,ti+1:end-1}};
              ei_detailed_data.TM(lsi, ti:end)         = ...
                {splitdd.TM{:}, ei_detailed_data.TM{lsi,ti+1:end-1}};
              ei_detailed_data.ei_info(lsi, ti:end)    = ...
                {splitdd.ei_info{:}, ei_detailed_data.ei_info{lsi,ti+1:end-1}};

              max_time_splits = max_time_splits + 1;
            end

          else
            error('Generation of collateral reduced space stopped after not converging in a time slice');
          end
        end

        % translate file-indices of snapshots into parameter vectors:
        ei_detailed_data.ei_info{lsi,ti-1}.extension_mus = ...
             Mtrain( : , ...
                     mod(ei_detailed_data.ei_info{lsi,ti-1}.extension_filenum-1, ...
                        length(Mtrain) )...
                       + 1 );
      end
    else % no adaptation
      params.ei_stop_epsilon = model.ei_stop_epsilon;
      for ti = 1:model.ei_time_splits
        params.time_index = ti;
        ei_detailed_data  = ei_detailed(model, ei_detailed_data, ...
                                        LU_fnames(lsi,:), params);
        % translate file-indices of snapshots into parameter vectors:
        ei_detailed_data.ei_info{lsi,ti}.extension_mus = ...
           Mtrain( : , ...
                   mod(ei_detailed_data.ei_info{lsi,ti}.extension_filenum-1, ...
                      length(Mtrain) )...
                     + 1 );
      end
    end

  end


 case 'file_load'
  disp('loading CRB from file for empirical interpolation');
  tmp = load(model.CRB_basis_filename);
  ei_detailed_data.QM = tmp.detailed_data.QM(:,1:model.Mmax);
  ei_detailed_data.BM = tmp.detailed_data.BM(:,1:model.Mmax);
  ei_detailed_data.TM = tmp.detailed_data.TM(1:model.Mmax);
  ei_detailed_data.ei_info = tmp.detailed_data.ei_info;
  if (model.Mmax~=size(tmp.detailed_data.QM,2))
    disp(['Only extracted ',num2str(model.Mmax),' CRB-vectors from'...
          num2str(size(tmp.detailed_data.QM,2))]);
  end;
  clear('tmp');

 otherwise
  error('ei-generation mode unknown');
end;

detailed_data.BM                 = ei_detailed_data.BM;
detailed_data.TM                 = ei_detailed_data.TM;
detailed_data.QM                 = ei_detailed_data.QM;
detailed_data.ei_info            = ei_detailed_data.ei_info;
detailed_data.time_split_map     = ei_detailed_data.time_split_map;

% test Mmax on consistency, otherwise there may be problems afterwards
if ~isequal(model.Mmax, size(detailed_data.BM{1},2))
  disp(['warning: Mmax in params is not the current number of generated', ...
        ' basis vectors. Please adjust model.Mmax!']);
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simply generate reduced basis as in lin-evol case

if model.Mstrich > 0
  model.M = cellfun(@(x)(size(x,2) - model.Mstrich), detailed_data.BM, 'UniformOutput', true);
end

detailed_data = rb_basis_generation(model,detailed_data);

