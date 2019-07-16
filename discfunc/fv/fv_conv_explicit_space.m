function INC = fv_conv_explicit_space(U,NU_ind,grid,params)
% function INC =
%      fv_conv_explicit_space(U,NU_ind,grid,params)
%  
% function applying an FV-space-discretization operator starting from old
% values U corresponding to the geometry given in grid producing a
% new vector of elementwise scalars NU but only on for the
% subelements with numbers given in NU_ind. If NU_ind is empty, all
% new values NU are determined, i.e. length(NU) = length(U) = grid.nelements
%
% By this, the operator evaluation can be performed in a localized
% way, i.e. used for empirical interpolation in rb_nonlin_evol_simulation
%
% usual timestepping can be performed afterwards by (NU = Id - deltat * INC). 
%
% required fields of params:
%            verbose:   a verbosity level

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


% Bernard Haasdonk 18.5.2007

% compute flux at dirichlet and inner edges. Neuman and
% cut-edges are not set.

%disp('step through fv_conv_explicit_space!!');
%keyboard;

num_flux = fv_num_conv_flux(U,grid,params);
num_flux_mat = num_flux.G;

% if no subset is specified: compute all elements
if isempty(NU_ind)
  NU_ind = 1:grid.nelements;
end;

%if params.verbose>=5
%  fprintf('.');
%end;

neu_NB_ind =find(grid.NBI == -2);
UU = repmat(U,1,grid.nneigh);

% determine neumann-boundary values at end as computed flux is required
% in case of outflow-conditions
if ~isempty(neu_NB_ind)
  
  % in case of file access, the correct filename must be set here
  % the following is only relevant in case of use of a 
  % params.use_velocitymatrix_file and filecaching mode 2
  if isfield(params,'filecache_velocity_matrixfile_extract') && ...
	(params.filecache_velocity_matrixfile_extract == 2);
    params.velocity_matrixfile = ... 
	cache_velocity_matrixfile_extract(...
	    grid.ECX(neu_NB_ind),...
	    grid.ECY(neu_NB_ind),...
	    'neumann_bnd', ...
	    params);
  end;
  
  FNneu = neuman_values( ...
      grid.ECX(neu_NB_ind),...
      grid.ECY(neu_NB_ind),...
      UU(neu_NB_ind),...
      grid.NX(neu_NB_ind),...
      grid.NY(neu_NB_ind), ...
      params);
  
  % set overall neumann boundary values
  num_flux_mat(neu_NB_ind) = grid.EL(neu_NB_ind) .* FNneu;  
end;

% for to_be_computed elements, we need all fluxes!
i = find(isnan(num_flux_mat(NU_ind,:)),1);
if ~isempty(i)
  error('not all fluxes specified, NaN occuring !!');
end;

%NU = U(NU_ind) - params.dt * grid.Ainv(NU_ind) .* ... 
%     sum( num_flux_mat(NU_ind), 2);  
INC = grid.Ainv(NU_ind) .* sum( num_flux_mat(NU_ind,:), 2);  

%  if params.verbose >= 10
%    plot_fv_data([U,NU],params,'U(n) and U(n+1)');
%    keyboard;
%  end;
% TO BE ADJUSTED TO NEW SYNTAX
%| \docupdate 
