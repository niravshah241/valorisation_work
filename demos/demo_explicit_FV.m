% demo of a explicit FV scheme for advection-diffusion
%
% Simple transport-diffusion problem: square domain with parabolic
% velocity profile, flow from left to right (or vice versa, depending on
% the sign of the specified velocity)
% upper half of the square is neuman-boundary with noflow conditions,
% lower half is dirichlet-boundary.
%
% A simulation for a square and a triangular grid is performed. For
% the triangular grid a more restrictive time-step is required, so
% 4 times as many time-slices are computed.

% Bernard Haasdonk 18.4.2006

help demo_explicit_FV

model = lin_evol_model_default;
gparams = {struct('gridtype','rectgrid'),...
           struct('gridtype','triagrid')};

%%% 1. common parameters
disp('Full FV-simulation:');
disp('-------------------');
disp('setting parameters');

model.xrange = [0,1]; % grid x-coordinate range
model.yrange = [0,1]; % grid y-coordinate range
model.xnumintervals = 100; % number of grid cells in x-direction
model.ynumintervals = 100; % number of grid cells in y-direction

% set everything to dirichlet (first column) 
% then the upper half boundary to neuman (second column)
model.bnd_rect_corner1 = [-1, -1; ...
                           -1,  0.5];
model.bnd_rect_corner2 = [ 2 , 2; ...
                            2 , 2];
model.bnd_rect_index =   [-1, -2]; % first is dirichlet, second neuman

% output options:
model.verbose = 20;  % 0   : no output, 
                      % >9  : informative output
                      % >19 : verbose informative output
                      % >29 : debug output 
                      % >99 : desperate debug output, with dbstops

model.T = 0.25;         % maximum time

%model.nt = 250;        % number of time steps

nt{1} = 125;% first for rectgrid, second for triagrid
nt{2} = 500;

% alternatively generate triagrid: smaller dt required for stability
%model.nt = model.nt*4;

% set data functions and simulation information
model.dirichlet_values_ptr = @dirichlet_values_homogeneous;
model.neumann_values_ptr   = @neumann_values_homogeneous;
model.c_dir = 1.0;
model.c_neu = 0.0;

% parameters inducing "affine dependency" into the problem
model.diffusivity_ptr   = @diffusivity_homogeneous;
model.num_diff_flux_ptr = @fv_num_diff_flux_gradient;
model.k = 0.004;      % diffusion coefficient 0 to 0.004 is OK
                       % although CFL condition is not completely
                       % met for highest value
model.laplacian_ptr            = @(glob, U, model) U;
model.laplacian_derivative_ptr = @(glob, U, model) ones(length(U),1);

model.conv_flux_ptr = @conv_flux_linear;
model.conv_flux_derivative_ptr = @conv_flux_forward_difference;
model.velocity_ptr  = @velocity_parabola;
model.flux_linear = 1;
model.num_conv_flux_ptr = @fv_num_conv_flux_engquist_osher;
model.c = -2;          % maximum velocity 0 to 2 is OK
model.flux_quad_degree = 1;

model.init_values_ptr = @init_values_blobs;
model.beta = 0.25;    % weighting of Gauss-blobs 0 to 1 is reasonable
model.gamma = 100;    % 100-1000 width of Gauss-blob of initial values
model.radius = 0.1;   % radius of cutoff-function

plot_params.show_colorbar = 1;
plot_params.no_lines      = 1;
plot_params.axis_equal    = 1;
plot_params.plot          = @plot_element_data;
model.mu_names = [];
model.debug = 10000;

%%% 2. "complete" simulation, complexities and plot

gridnames = {'rectangular','triangular'};
fluxes = {@fv_num_conv_flux_lax_friedrichs,...
          @fv_num_conv_flux_engquist_osher};
for g = 1:2
  disp(['starting full simulation for ',gridnames{g},' grid and ',...
        ' numerical flux ',func2str(fluxes{g})]);
  tic;
  model.num_conv_flux_ptr = fluxes{g};
  model.nt                = nt{g};
  model                   = model_default(model);
  model                   = structcpy(model, gparams{g});
  model_data = gen_model_data(model);
  if g == 2
    % use standard unit-square grid ...
    model_data.grid = triagrid;
    % ... and update boundary data
    model_data.grid = set_boundary_types(model_data.grid, model);
  end
  U = fv_conv_diff(model,model_data);
  t = toc;
  disp(['===> Computation time of full simulation : ',num2str(t)]);
  disp(['===> Memory requirement for full simulation : ', ...
       num2str(numel(U)*8/(1024*1024)), ' MB']);
  % U = rand(model.nx*model.ny,model.nt+1);

  model.title = ['Simulation ',gridnames{g},' grid, ',...
                 func2str(fluxes{g}),' numerical flux.'];

  plot_sequence(U, model_data.grid, plot_params);
end;


