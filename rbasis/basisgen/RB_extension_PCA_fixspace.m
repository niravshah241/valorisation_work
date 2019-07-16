function [RBext, dummy] =  RB_extension_PCA_fixspace(model, ...
                                                     detailed_data)
%function [RBext, dummy] = RB_extension_PCA_fixspace(model, ...
%                                                   detailed_data)
% function computing a RB basis extension for given parameters by
% the POD-Greedy algorithm.
%
% The choice is made by complete time simulation for the current parameter `mu`
% set via last call to set_mu() and taking the principal component of PCA
% keeping the RB as fixspace. This single vector is orthonormalized with
% respect to RB and returned in RBext.
%
% Return argument dummy is a superfluous argument not used in this routine, but
% necessary for a uniform argument list compatible with other extension
% algorithms.
%
% required fields of model:
%  mu_names                 : cell array of names of the parameter-vector
%                             entries
%
% return values:
%  RBext : the new reduced basis vector meant for extension of the current
%          reduced basis space
%  dummy : dummy variable making argument list compatible with other extenstion
%          algorithms.
%

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

% Bernard Haasdonk 14.10.2005

dummy = [];

% cache detailed simulation
%  sim_data = cached_detailed_simulation(model,detailed_data,sim_params);
%model_data.grid = detailed_data.grid;
%model_data.W    = detailed_data.W;

%tmp_params.mu = sim_params.mu;
%model = model.set_mu(model,mu);
% use detailed_data as model_data here

% take care, that the arguments do not change throughout parameter variation!
% detailed_data: RB_info, V, W
% model: N, Nmax, base_model


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

sim_data = filecache_function(@detailed_simulation,tmp_model, ...
                              tmp_detailed_data);

U = model.get_dofs_from_sim_data(sim_data);

% select single vector for extension
%keyboard;
RBextall = PCA_fixspace(U,...
    model.get_rb_from_detailed_data(detailed_data), ...
    model.get_inner_product_matrix(detailed_data), ...
    1);

if ~isempty(RBextall)
  RBext = RBextall(:,1);
else
  disp('warning: did not find basis extension in PCA_fixspace!!!');
  RBext = [];
end;

%disp('check new basis-vector!!')
%keyboard;

