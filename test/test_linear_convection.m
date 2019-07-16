function OK = test_linear_convection
% function simulating a simple linear convection scheme with source term
% u_t + b*Du + cu = 0,               u(x,0) = u0.
% u0 is 1 inside a box and 0 otherwise. This box gets transported in
% direction b with a slowing term c.
% returns 1, if error to exact solution is "small enough".
% returns 0 otherwise

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


% Martin Drohmann 26.05.2008

OK = 1;

params = [];

% rectgrid or trigrid
params.gridtype = 'rectgrid';
% params.grid_initfile = '2dsquaretriang';

%params.gridtype = 'rectgrid';
params.xnumintervals = 40;
params.ynumintervals = 40;
params.xrange = [0,1];
params.yrange = [0,1];

params.bnd_rect_corner1 = [-1, -1];
params.bnd_rect_corner2 = [2, 2];

params.bnd_rect_index = (-2);
params.name_neuman_values = 'outflow';

params.T = 1.0;
params.nt = 400;

params.axis_equal = 1;

params.name_init_values = 'gradient_box';
%params.implict_operators_algorithm = 'eye';
params.c_init_up = 1.0;
params.c_init_lo = 1.0;

params.name_flux = 'transport';
params.flux_linear = 1;
params.verbose = 1;
params.debug = false;
params.transport_x = 0.1;
params.transport_y = 0.1;
params.transport_source = 1;
params.no_lines = true;

params.name_convective_num_flux  = 'lax-friedrichs';
params.convective_discretization = 'explicit';
params.diffusive_discretization  = 'none';
params.reaction_discretization   = 'explicit';

% for a detailed simulation we do not need the affine parameter
% decomposition
params.affine_decomp_mode = 'complete';

grid = construct_grid(params);

% initial values by midpoint evaluation
U0 = fv_init_values(grid,params);

U = zeros(length(U0(:)),params.nt+1);
U(:,1) = U0(:);
params.dt = params.T/params.nt;
dt = params.dt;

params.data_const_in_time = 1;

% loop over fv-steps
for t = 1:params.nt  
  params.t = t;
  if params.verbose > 19
    disp(['entered time-loop step ',num2str(t), ' for transport equation']);
  elseif params.verbose > 9
    fprintf('.');
  end;

  [inc, fluxes] = fv_explicit_space(U(:,t), [], ...
              grid, params);

  U(:,t+1) = U(:,t) - params.dt * inc;
end

EX = zeros(length(U0(:)), params.nt + 1);
EX(:,1) = U0(:);

for t = 1:params.nt
  EX(:,t+1) = exp(-params.transport_source * t * dt ) * ...
      init_values(grid.CX(:)-params.transport_x * t *dt, ...
      grid.CY(:) - params.transport_y * t * dt, params);
end

error = fv_l2_error(U(:,params.nt+1), EX(:,params.nt+1), grid, params);
if error > 0.1
    OK = 0;
end
% plot_element_data_sequence(grid,U,params);
% plot_element_data_sequence(grid,EX,params);

% TO BE ADJUSTED TO NEW SYNTAX
%| \docupdate 
