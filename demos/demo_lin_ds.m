function demo_lin_ds
%function demo_lin_ds
%
% function demonstrating RB approach for linear dynamical systems
% according to MATHMOD 2009 paper
%
% Method also demonstrates new call of all reduced schemes:
%
% model_data = gen_model_data(model);
%
% sim_data = detailed_simulation(model,model_data); 
%    produces some structure with suitable fields, e.g DOF vector
%    output sequence, reference to output files, simulation
%    information, etc. Last parameter is model-struct, possibly
%    precedent parameters, e.g. grid, etc.
%    so not only DOF-vector is the return argument, but a structure
% detailed_data     = gen_detailed_data(model,model_data); 
%    produce high-dimensional data, e.g. reduced basis, etc.
%    mu is not known.
% reduced_data      = gen_reduced_data(model,detailed_data); 
%    produce low-dimensional data, e.g. gram matrices between
%    reduced basis vectors, subgrid extraction, etc.
%    mu is not known (so is the previous offline-data or
%    combination of offline and online data)
% rb_sim_data  = rb_simulation(model,reduced_data)
%    perform reduced simulation, i.e. mu known
%    produces some structure with suitable fields, e.g DOF vector,
%    error estimation
%
% Further routines "extending" datasets, e.g.:
%     rb_sim_data   = rb_reconstruction(rb_sim_data,model)
%     reduced_data   = extract_reduced_data_subset(reduced_data,model)
%     detailed_sim_data   = rb_basis_generation(detailed_data,model)

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


% Bernard Haasdonk 2.4.2009

%step = 1 % simulation and visualization of dynamical system
%step = 2 % reduced simulation and visualization of approximation
step =3 % detailed and reduced simulation and visualization of
        % error and estimator

	
model = lin_ds_model_default;
	
model.data_const_in_time = 1;
model.verbose = 0;
model.T = 5*pi;
model.nt = 10000;
model.affinely_decomposed = 1;
model.A_function_ptr = @A_func;
% ratio of first to second ellipses axes
%model.A_axes_ratio = 2;
model.A_axes_ratio = 2;
model.B_function_ptr = @B_func;
model.C_function_ptr = @C_func;
model.D_function_ptr = @D_func;
model.x0_function_ptr = @x0_func;
% shift of x0 in third coordinate
model.x0_shift=1;
model.u_function_ptr = @(model) 0.1; % define u constant 
%model.u_function_ptr = @(model) sin(4*model.t);  
%model.u_function_ptr = @(model) 0; % define u constant 0
%model.u_function_ptr = @(model) model.t; % define u = 1 * t
model.G_matrix_function_ptr = @(model,model_data) eye(3,3);
% the following lead to both larger C1 and C2!!!
%model.G_matrix_function_ptr = @(model) diag([2,0.5,1]);
%model.G_matrix_function_ptr = @(model) diag([0.5,2,1]);

% dummy parametrization:, parameter is not used...
model.mu_names = {'A_axes_ratio','x0_shift'};
%model.mu_ranges = {[0.5,2],[0,1]};
model.mu_ranges = {[0.5,0.5],[1,1]};
model.affinely_decomposed = 1;
 
% set function pointers to main model methods
%model.gen_model_data = @lin_ds_gen_model_data;
%model.detailed_simulation = @lin_ds_detailed_simulation;
%model.gen_detailed_data = @lin_ds_gen_detailed_data;
%model.gen_reduced_data = @lin_ds_gen_reduced_data;
%model.rb_simulation = @lin_ds_rb_simulation;
%model.plot_sim_data = @lin_ds_plot_sim_data;
%model.rb_reconstruction = @lin_ds_rb_reconstruction;
% determin model data, i.e. matrix components, matrix G

model_data = gen_model_data(model); % matrix components

% be sure to determine correct constants once, if anything in the
% problem setting has changed!!
estimate_bounds = 1;
% the following only needs to be performed, if the model are changed.
if estimate_bounds
  model.estimate_lin_ds_nmu = 10;
  model.estimate_lin_ds_nX = 100;
  model.estimate_lin_ds_nt = 100;
  format long;
  [C1,C2] = lin_ds_estimate_bound_constants(model,model_data)
  model.state_bound_constant_C1 = C1;
  model.output_bound_constant_C2 = C2;
else
  model.state_bound_constant_C1 = 2;
  model.output_bound_constant_C2 = sqrt(3);
end;
model.error_estimation = 1;
model_data.RB = eye(3,2); 
model.RB_generation_mode = 'from_detailed_data';

switch step
 case 1 % simulation and visualization
  %  later: sim_data = detailed_simulation(model);
    
  sim_data = detailed_simulation(model,model_data);

  p = plot_sim_data(model,model_data,sim_data,[]);
  
 case 2 % reduced simulation
  %define reduced basis of first two coordinates (i.e. error only
  %by u!)

  model.debug = 1;
  
  %    produce high-dimensional data, e.g. reduced basis, etc.
  %    mu is not known
  detailed_data     = gen_detailed_data(model,model_data); 
  
  %    produce low-dimensional data, e.g. gram matrices between
  %    reduced basis vectors, subgrid extraction, etc.
  %    mu is not known (should comprise previous offline and online steps)

  reduced_data = gen_reduced_data(model,detailed_data);
  model.N = reduced_data.N;
  
  % rb_sim_data  = rb_simulation(reduced_data,model)
  %    perform reduced simulation, i.e. mu known
  %    produces some structure with suitable fields, e.g DOF vector
  
  rb_sim_data = rb_simulation(model,reduced_data);
  
  rb_sim_data = rb_reconstruction(model,detailed_data,rb_sim_data);

  model.plot_title = 'reduced trajectory';
  %model.plot_indices = 1:3;   => default
  p = plot_sim_data(model,model_data,rb_sim_data,[]);
%  keyboard;

 case 3 % detailed and reduced simulation
  model.x0_shift=0;
  %model.u_function_ptr = @(model) 0; % define u constant 0
  % => leads to error 0 :-) if x0_shift = 0
  
  %model.u_function_ptr = @(model) 0.1; % define u constant 
  %=> error linearly increasing
  model.u_function_ptr = @(model) 0.3*sin(2*model.t)+0.002*model.t;     

  
  sim_data = detailed_simulation(model,model_data);
  detailed_data     = gen_detailed_data(model,model_data); 
  model.enable_error_estimator=1;
  reduced_data = gen_reduced_data(model,detailed_data);
  model.N = reduced_data.N;
  rb_sim_data = rb_simulation(model,reduced_data);
  rb_sim_data = rb_reconstruction(model,detailed_data, ...
					 rb_sim_data);
  
  figure, plot3(sim_data.X(1,:),...
		sim_data.X(2,:),...
		sim_data.X(3,:), ...
                'color','b' ...
		);
  hold on
  plot3(rb_sim_data.X(1,:),...
		rb_sim_data.X(2,:),...
		rb_sim_data.X(3,:),...
		'color','r');
  axis equal;
  title('exact and reduced trajectories');
  legend({'exact','reduced'});
  
  error = sim_data;
  error.X = sim_data.X-rb_sim_data.X;
  error.Y = sim_data.Y-rb_sim_data.Y;

  err = sqrt(sum((detailed_data.G * error.X).*error.X,1)); 
  figure, plot([rb_sim_data.DeltaX(:),err(:)]), 
  Ylim = get(gca,'Ylim');
  set(gca,'Ylim',[0,max(Ylim)]);
  title('rb state error and estimator');
  legend({'estimator','error'});

  figure, plot([rb_sim_data.DeltaY(:),abs(error.Y(:))]); 
  title('rb output error and estimator');
  Ylim = get(gca,'Ylim');
  set(gca,'Ylim',[0,max(Ylim)]);
  legend({'estimator','error'});
  disp('end');
  %keyboard;
  
 otherwise
  error('step unknown');
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% functions defining the dynamical system and input/initial data
function A = A_func(model,model_data);
A = eval_affine_decomp(@A_components,@A_coefficients,...
		       model,model_data);

% A_parameter: ratio of first to second axis of 
function Acomp = A_components(model,model_data)
%Acomp = {[0 -0.25 0; 4 0 0 ; 0 0 0 ]};
%Acomp = {[0 -1 0; 1 0 0 ; 0 0 0 ]};
%Acomp = {[-0.01 0 0; 0 -0.9 0 ; 0 0 -0.001 ]};
Acomp = {[0 1 0; 0 0 0 ; 0 0 0 ],...
	 [0 0 0; 1 0 0 ; 0 0 0 ] ...
	};

function Acoeff = A_coefficients(model)
Acoeff = [-model.A_axes_ratio, 1/model.A_axes_ratio];

function B = B_func(model,model_data)
B = eval_affine_decomp(@B_components,@B_coefficients,model,model_data);

function Bcomp = B_components(model,model_data)
Bcomp = {[0;0;1]};

function Bcoeff = B_coefficients(model)
Bcoeff = 1;

function C = C_func(model,model_data)
C = eval_affine_decomp(@C_components,@C_coefficients,model,model_data);

function Ccomp = C_components(model,model_data)
Ccomp = {[1,1,1]};
%Ccomp = {[2,0,1],[0,2,1]};

function Ccoeff = C_coefficients(model)
Ccoeff = 1;
%Ccoeff = [model.C_weight,1-model.C_weight];

function D = D_func(model,model_data)
D = 0; % no decomposition required for scheme

function x0 = x0_func(model,model_data)
x0 = eval_affine_decomp(@x0_components,@x0_coefficients,model,model_data);

function x0comp = x0_components(model,model_data)
x0comp = {[1;0;0],[1;0;1]};

function x0coeff = x0_coefficients(model)
x0coeff = [1-model.x0_shift,model.x0_shift];




%| \docupdate 
