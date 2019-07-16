function [RBext,info] =  RB_extension_max_error_snapshot(model,detailed_data)
%function [RBext,par] =  RB_extension_max_error_snapshot(model,detailed_data)
%
% function computing the snapshot for parameter mu which is set in
% params, which causes the largest increase of the posterior error estimator 
% in time. This snapshot is orthonormalized with respect to detailed_data.RB and
% returned in RBext.
%
% Occasionally, the selected snapshot is already within the span
% of the reduced basis RB, then the second max-error snapshot is taken, 
% and so on.
%
% The values returned in info are:
%      snapshot: the snapshot before orthonormalization
%      ti: the absolute simulation time of the snapshot
%      rank: rank of the snapshot in the list of 
%            error-prediction.increases
%
% required fields of params 
%      detailed_simulation_algorithm: name of the function to be called
%           with parameters (grid,params) and resulting in DOF-vectors
%           e.g. fv_conv_diff_operators
%           or detailed_simulation with additional fields
%      mu_names:   cell array of names of the parameter-vector entries
%      inner_product_matrix_algorithm : function computing the correct inner
%                product matrix between function DOF-vectors

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


% Bernard Haasdonk 13.10.2005

  info = []; 

  % perform full and RB simulation
  %U = fv_conv_diff(params);
  %if isempty(grid)
  %  grid = grid_geometry(params);
  %end;

  % take care to remove some changed fields from model and detailed
  % data, such that filecaching works!!! in RB_extension_algorithm
  % !!!
  tmp_model = model;
  ifm = model.filecache_ignore_fields_in_model;
  for i = 1:length(model.filecache_ignore_fields_in_model);
    if isfield(model,ifm{i})
      tmp_model = rmfield(tmp_model,ifm{i});
    end;
  end;
  ifm = model.filecache_ignore_fields_in_detailed_data;
  dim_RB = size(model.get_rb_from_detailed_data(detailed_data));
  tmp_detailed_data = ...
      model.set_rb_in_detailed_data(detailed_data, zeros(dim_RB(1),0));
  for i = 1:length(model.filecache_ignore_fields_in_detailed_data);
    if isfield(detailed_data,ifm{i})
      tmp_detailed_data = rmfield(tmp_detailed_data,ifm{i});
    end;
  end;
  %  U = fv_conv_diff_operators(grid,params);
  sim_data = filecache_function(@detailed_simulation,tmp_model, ...
                                tmp_detailed_data);

  % perform RB-simulation
  reduced_data = gen_reduced_data(model,detailed_data);
  simulation_data = rb_simulation(model,reduced_data);
  Urec         = rb_reconstruction(model,detailed_data,simulation_data);
  %[A,pe] = fv_conv_diff_galerkin(RB,params);
  %Uappr = RB*A;
  simulation_data.Delta = fv_error(Urec,U,detailed_data.grid,model);

  % rank timesteps according to error increase
  post_est_diff = simulation_data.Delta' - ...
      [0; simulation_data.Delta(1:end-1)'];
  [m,ti] = sort(-post_est_diff); % -> post_est_diff(ti) ist sortiert

  % select snapshot by following rank-list
  ra = 0;
  found_update = 0;

  W = model.get_inner_product_matrix(detailed_data);

  while ~found_update
    ra = ra+1;
    % orthonormalize snapshot at rank ra
    Uproj = U(:,ti(ra)) - detailed_data.RB * (detailed_data.RB' * W * U(:,ti(ra)));
    if norm(Uproj)> 1e-10 
      % if nonzero projection error -> end with new snapshot
      %    Uproj = Uproj/norm(Uproj).*params.vec_to_func;
      % due to scaling, no precise orthonormalization is obtained. So perform
      % explicit orthonormalization
      RBtmp = orthonormalize([detailed_data.RB,Uproj],W);
      if size(RBtemp,2)> size(detailed_data.RB,2)
        Uproj = RBtmp(:,end); %.*grid.Ainv(:).^(0.5);
        found_update = 1;
      else
        error('accuracy problem, should not reach this point!!');
      end;
    end;
  end;

  % return values
  RBext = Uproj;
  info.rank = ra;
  info.snapshot = U(:,ra);
  info.ti = (ti(ra)-1)*params.T/params.nt;

%| \docupdate 
