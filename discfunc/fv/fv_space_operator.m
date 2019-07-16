function [INC, fluxes] = fv_space_operator(model,model_data,U,NU_ind,weights)
% function INC =
%      fv_space_operator(model,model_data,U,NU_ind,grid,weights)
%
% function applying an FV-space-discretization operator starting from old
% values U corresponding to the geometry given in model_data.grid producing a
% new vector of elementwise scalars NU but only on for the subelements with
% numbers given in NU_ind. If NU_ind is empty, all new values NU are
% determined, i.e. length(NU) = length(U) = grid.nelements
%
% By this, the operator evaluation can be performed in a localized
% way, i.e. used for empirical interpolation in e.g. rb_nonlin_evol_simulation
%
% usual timestepping can be performed afterwards by (NU = Id - deltat *
% INC).
%
% required fields of model:
%            verbose:   a verbosity level
% model shall set the fields:
%            diff_weight:    weight for diffusive flux component
%            conv_weight:    weight for convective flux component
%            reac_weight:    weight for reaction term

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


% Martin Drohmann 9.12.2007 based on fv_conv_explicit_space by Bernard
% Haasdonk

% compute flux at dirichlet and inner edges. Neuman and
% cut-edges are not set.

grid = model_data.grid;

if ~model.debug
  warning('off', 'MATLAB:structOnObject');
end
smodel = struct(model);

num_diff_flux.G=zeros(size(grid.NBI));
num_conv_flux.G=zeros(size(grid.NBI));
if weights.diff_weight ~= 0
  num_diff_flux = model.num_diff_flux_ptr(smodel,model_data,U,NU_ind);
  num_flux_mat  = weights.diff_weight * num_diff_flux.G;
else
  num_flux_mat = zeros(size(grid.NBI));
end
if weights.conv_weight ~= 0
  num_conv_flux = model.num_conv_flux_ptr(smodel,model_data,U);
  num_flux_mat  = num_flux_mat + weights.conv_weight * num_conv_flux.G;
end

% if no subset is specified: compute all elements
if isempty(NU_ind)
  NU_ind = 1:grid.nelements;
end

%if model.verbose>=5
%  fprintf('.');
%end;

neu_NB_ind = find(grid.NBI == -2);
UU = repmat(U,1,grid.nneigh);

% determine neumann-boundary values at end as computed flux is required
% in case of outflow-conditions
if ~isempty(neu_NB_ind>0)
  % in case of file access, the correct filename must be set here
  % the following is only relevant in the case where a
  % model.use_velocitymatrix_file and filecaching mode 2 are selected
  if model.filecache_velocity_matrixfile_extract == 2;
        model.velocity_matrixfile = ... 
          velocity_matrixfile_extract(...
              model, ...
              'neumann_bnd', ...
	      grid.ECX(neu_NB_ind),...
	      grid.ECY(neu_NB_ind));
  end

  FNneu = model.neumann_values_ptr( ...
      [grid.ECX(neu_NB_ind), grid.ECY(neu_NB_ind)],...
      UU(neu_NB_ind),...
      [grid.NX(neu_NB_ind),grid.NY(neu_NB_ind)], model);
  
  % set overall neumann boundary values
  num_flux_mat(neu_NB_ind) = grid.EL(neu_NB_ind) .* FNneu;
end;

% for to_be_computed elements, we need all fluxes!
if ~isempty(find(isnan(num_flux_mat(NU_ind,:)),1))
    disp(['mu = ( ', mat2str(get_mu(model)), ' ) ']);
  error('not all fluxes specified, NaN occuring !!');
end;

%NU = U(NU_ind) - model.dt * grid.Ainv(NU_ind) .* ... 
%     sum( num_flux_mat(NU_ind), 2); 
INC = grid.Ainv(NU_ind,:) .* sum( num_flux_mat(NU_ind,:), 2);

react = zeros(1,grid.nelements);
if weights.react_weight ~= 0
  react = reaction(model, grid.CX, grid.CY, U);
  INC   = INC + weights.react_weight * react(NU_ind,:);
end
fluxes = [];

%  if model.verbose >= 10
%    plot_fv_data([U,NU],model,'U(n) and U(n+1)');
%    keyboard;
%  end;
