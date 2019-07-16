function NU = fv_conv_diff_step(U,num_conv_flux,num_diff_flux,grid,params)
% function NU = fv_conv_diff_step(U,num_conv_flux,num_diff_flux,grid,params)
%
% function performing an FV-time-step starting from old
% values U corresponding to the geometry given in grid producing a
% new vector of elementwise scalars NU using given flux-information
% neuman-flux treatment is performed, dirichlet is assumed to be
% incorporated in the specified fluxes num_conv_flux, num_diff_flux.

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


% required fields of params:
% dt        : size of the timestep
% k         : diffusion parameter

% Bernard Haasdonk 27.9.2005

  num_flux_mat = num_conv_flux.G + num_diff_flux.G;

  if params.verbose>=5
    fprintf('.');
  end;

  neu_NB_ind =find(grid.NBI == -2);
  UU = repmat(U,1,size(grid.NBI,2));

  % determine neumann-boundary values at end as computed flux is required
  % in case of outflow-conditions
  if ~isempty(neu_NB_ind)
    FNneu = params.neumann_values_ptr( ...
        [grid.ECX(neu_NB_ind),...
         grid.ECY(neu_NB_ind)],...
        UU(neu_NB_ind),...
        [grid.NX(neu_NB_ind),...
         grid.NY(neu_NB_ind)], ...
        params);

    % set overall neumann boundary values
    num_flux_mat(neu_NB_ind) = grid.EL(neu_NB_ind) .* FNneu;
  end;

  i = find(isnan(num_flux_mat),1);
  if ~isempty(i)
    error('not all fluxes specified, NaN occuring !!');
  end;

  NU = U - params.dt * grid.Ainv .* sum( num_flux_mat , 2);

%  if params.verbose >= 10
%    plot_fv_data([U,NU],params,'U(n) and U(n+1)');
%    keyboard;
%  end;
%| \docupdate 
