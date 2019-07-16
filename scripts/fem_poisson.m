function res = fem_poisson(step,params)
%function res = fem_poisson(step)
%
% Script demonstrating lagrange finite element functions, using
% them for function interpolation and finite element methods.
%
% Script realizing the simple poisson equation or a more complex
% elliptic problem with arbitrary diffusivity_tensor, advection,
% source, reaction, neumann , dirichlet and robin boundary
% conditions on the unit square
% with finite elements on triagrid discretization.
%
% step 1: some basic femdiscfunc demonstration and functionatity
% step 2: use femdiscfunc for interpolating functions by local and
%         global evaluations
% step 3: solve fem-problem, return errors
% step 4  error convergence investigation
% step 5: same as 3 by call of interface methods + Reduced basis steps
% step 6: demo_rb_gui

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


% B. Haasdonk 11.1.2011

% adapted to new assembly

% Immanuel Maier 12.04.2011

% poisson_model

if nargin == 0
  step = 1;
end;

if nargin < 2
  params = [];
end;

res = [];

switch step
  
 case 1 % elementary fem discrete function operations

  params = [];
  pdeg = 4;
  params.pdeg = pdeg;
  params.dimrange = 3;
  % params.dimrange = 1;
  params.debug = 1;
  params.numintervals = 5;
  model = poisson_model(params);
  % convert to local_model:
  model = elliptic_discrete_model(model);
  grid = construct_grid(model);
  df_info = feminfo(params,grid);
  %tmp = load('circle_grid');
  %grid = triagrid(tmp.p,tmp.t,[]);
  disp('model initialized');

  % without model_data, detailed_simulation, etc. but explicit
  % calls of fem operations
    
  % initialize vectorial discrete function, extract scalar
  % component and plot basis function
  df = femdiscfunc([],df_info); % initialize zero df
  
  fprintf('\n');
  disp('initialization of femdiscfunc successful, display:');
  display(df);
  
  disp('press enter to continue');
  pause;
  
  % for check of dof consistency we have a plot routine:
  fprintf('\n');
  disp('plot of global_dof_index map for consistency check:');
  p = plot_dofmap(df);
  
  disp('press enter to continue');
  pause;
  
  % global dof access:
  fprintf('\n');
  disp('example of dof access, scalar component extraction and plot:');    
  % set second vector component in 4th basis function nonzero;
  ncomp = 2;
  locbasisfunc_index = 4;
  df.dofs((locbasisfunc_index-1)*df.dimrange+ncomp) = 1; 
  dfscalar = scalar_component(df,2); % should be nonzero function
  disp(['entry should be 1 : ',...
	num2str(dfscalar.dofs(locbasisfunc_index))]); 
  
  params = [];
  params.subsampling_level = 6;
  figure,plot(dfscalar,params);
  dfscalar.dofs(:) = 0; % now zero again.
%  keyboard;
  
  disp('press enter to continue');
  pause;
  
  fprintf('\n');
  disp('example of local dof access:');
  % local dof access via local to global map
  eind = 4; % set dof on element number 4
  lind = (pdeg+2); % local basis function index with lcoord (0,1/pdeg) 
  dfscalar.dofs(dfscalar.global_dof_index(eind,lind)) = 1;
  
  % example of local evaluation of femdiscfunc simultaneous on
  % several elements in the same local coordinate point
  elids = [4,6,7,10]; % evaluation on some elements
  lcoord = [0,1/pdeg]; % local coordinate vector == last point in all triangles
  f = evaluate(dfscalar,elids,lcoord); 
  disp(['first entry should be 1 : ',num2str(f(1))]); 
  % equivalent call (!) by () operator as abbreviation for local evaluation:
  f = dfscalar(elids,lcoord); 
  disp(['first entry should be 1 : ',num2str(f(1))]); 
  
  disp('press enter to continue');
  pause;

  disp('examples of norm computation:')
  params.dimrange = 1;
  params.pdeg = 1;
  dfinfo1 = feminfo(params,grid);
  df1 = femdiscfunc([],dfinfo1);
  df1.dofs(:) = 1;
  disp(['L2-norm(f(x,y)=1) = ',num2str(fem_l2_norm(df1))]);
  disp(['H10-norm(f(x,y)=1) = ',num2str(fem_h10_norm(df1))]);
  df1.dofs(:) = df1.grid.X(:);
  disp(['L2-norm(f(x,y)=x) = ',num2str(fem_l2_norm(df1))]);
  disp(['H10-norm(f(x,y)=x) = ',num2str(fem_h10_norm(df1))]);
  df1.dofs(:) = df1.grid.Y(:);
  disp(['L2-norm(f(x,y)=y) = ',num2str(fem_l2_norm(df1))]);
  disp(['H10-norm(f(x,y)=y) = ',num2str(fem_h10_norm(df1))]);
  
  disp('press enter to continue');
  pause;
  
  % evaluate df in all lagrange nodes of element 4 by loop
  fprintf('\n');
  disp(['dfscalar on element 4 in all lagrange nodes,' ...
	'only (pdeg+2) entry should be 1:']);
  lagrange_nodes = lagrange_nodes_lcoord(pdeg);
  elid = 4;
  for i = 1:size(lagrange_nodes,1);
    f = evaluate(dfscalar,elid,lagrange_nodes(i,:)); 
    disp(['f(l(',num2str(i),')) = ',num2str(f)]);
  end;

  disp('press enter to continue');
  pause;
  
  fprintf('\n');
  disp('example of requirement of subsampling in plot of discfuncs:');
  figure;
  subsamp_levels = [0,2,4,16];
  for i=1:length(subsamp_levels)
    subplot(2,2,i),
    params.axis_equal = 1;
    params.subsampling_level = subsamp_levels(i);
    params.clim = [-0.15    1.15]; % rough bounds
    plot(dfscalar,params);
    title(['subsampling level = ',num2str(subsamp_levels(i))]);
  end;

  disp('end of step 1, please inspect variables, if required.')
  keyboard;    

 case 2 

  params = [];
  pdeg = 4;
  params.pdeg = pdeg;
  params.dimrange = 1;
  params.numintervals = 5;
  model = poisson_model(params);
  % convert to local_model:
  model = elliptic_discrete_model(model);
  grid = construct_grid(model);
  %tmp = load('circle_grid');
  %grid = triagrid(tmp.p,tmp.t,[]);
  disp('model initialized');
  
  % interpolate exact solution and other coefficient functions 
  % as fem-function and plot

  disp('examples of interpolation of analytical functions:');
  df_info=feminfo(model,grid);
  df = femdiscfunc([],df_info);

  df = fem_interpol_global(model.solution, df);  
  plot(df);
  title('analytical solution');

  % discretize source function and plot

  % problems with fem_interpol_local!!
  %df = fem_interpol_local(model.source, df);
  df = fem_interpol_global(model.source, df); 
  figure,plot(df);
  title('source function');

  disp('press enter to continue');
  pause;
  
  % discretize diffusivity and plot as 4-sequence
  % to be implemented for vectorial functions... 
  %params.dimrange = 4;
  %df4 = femdiscfunc([],grid,params);
  %df4 = fem_interpol_global(model.diffusivity_tensor, df)  
  % arrange as sequence of 4 scalar functions and plot
  %plot(df4);
  %title = 'diffusivity function';
  
  % example of local evaluation of vectorial basis function 
  % on reference triangle
  disp('local evaluation of vectorial basis functions on reference triangle:');
  lcoord = [0,1];
  params = [];
  params.dimrange = 3;
  params.pdeg = 2;
  df2 = femdiscfunc([],df_info);
  res = fem_evaluate_basis(df2,lcoord)

  disp('local evaluation of scalar basis functions derivative on reference triangle:');
  res = fem_evaluate_scalar_basis_derivative(df,lcoord)
  
  % later:
  % disp('local evaluation of vectorial basis function derivative on reference triangle:');

  % later: 
   % interpolation error convergence with grid refinement
 
 case 3 % finite element assembly, solution and plot.
  % different components assembled separately, as later
  % affine-decomposition will require such modularization
  
  % assemble discrete system, solve and plot results

  if isempty(params)
    params = [];
    params.dimrange = 1;
    pdeg = 1;
    qdeg = 3 * pdeg; % quadrature_degree
    params.pdeg = pdeg;
    params.qdeg = qdeg;
    params.numintervals = 20;
    params.no_plot = 0; %plotting active
  end;

  no_plot = params.no_plot;
  
  %model = poisson_model(params);
  model = elliptic_debug_model(params);
  model = elliptic_discrete_model(model); % convert to local_model
  grid = construct_grid(model);
  disp(['grid.nelements=',num2str(grid.nelements)]);
  % for debugging: set all to neumann
  %  i = find(grid.NBI==-1);
  %  grid.NBI(i)=-2;
  
  df_info = feminfo(model,grid);
  disp(['df_info.ndofs=',num2str(df_info.ndofs)]);
  %tmp = load('circle_grid');
  %grid = triagrid(tmp.p,tmp.t,[]);
  disp('model initialized');
  
  %%%%%%%%%%%%%%% assembly of system %%%%%%%%%%%%%%%%%%%

  % assemble right hand side
  [r_source, r_dirichlet, r_neumann, r_robin] = ...
      fem_rhs_parts_assembly(model,df_info);
  r = r_source + r_neumann + r_robin + r_dirichlet;
  disp('rhs assembled');
  
  % sparse system matrix:
  A = spalloc(df_info.ndofs,df_info.ndofs,10); % rough upper bound for number of nonzeros
  [A_diff , A_adv, A_reac, A_dirichlet, A_neumann, A_robin] = ...
      fem_matrix_parts_assembly(model,df_info);  
  A = A_diff + A_adv + A_reac + A_neumann + A_robin + A_dirichlet;
  disp('matrix assembled');
    
  %%%%%%%%%%%%%%% solve system %%%%%%%%%%%%%%%%%%%
  % solution variable
  u_h = femdiscfunc([],df_info);
  u_h.dofs = A\r;    
  disp('system solved');
    
  %%%%%%%%%%%%%%% postprocessing %%%%%%%%%%%%%%%%%%%

  if ~no_plot
    disp('plotting...');
    % plot result
    params.subsampling_level = 10;
    figure, plot(u_h,params);
    title('discrete solution')
  end;   
  
  df = femdiscfunc([],df_info);
  df = fem_interpol_global(model.solution, df,model);  
  if ~no_plot
    figure, plot(df);
    title('analytical solution');
  end;
  
  % plot error of (interpolated) exact solution and discrete sol.
  % interpolation using same degree as discrete solution
  e_h = femdiscfunc([],df_info);
  e_h = fem_interpol_global(model.solution, e_h,model);  
  e_h.dofs = e_h.dofs -u_h.dofs;
  if ~no_plot
    figure, plot(e_h);
    title('rough error I\_rough(u)-u_h');
  end;
  
  if ~no_plot
    % better: choose high degree for interpolation of u
    % then require local_interpolation of u_h into high-resolved space
    params.pdeg = 4; % choose high degree
		     
    % map u_h to higher degree discretefunction
    params.has_dirichlet_values = 10;
    df_info_fine = feminfo(params,grid);
    u_h_fine = femdiscfunc([],df_info_fine);
    u_h_fine = fem_interpol_local(u_h, u_h_fine,model);  
    e_h_fine = femdiscfunc([],df_info_fine);
    e_h_fine = fem_interpol_global(model.solution, e_h_fine,model);  
    e_h_fine.dofs = e_h_fine.dofs -u_h_fine.dofs;
    figure, plot(e_h_fine);
    title('fine error I\_fine(u)-u_h');
  end;
  
  % determine error in L2 and H10 norm
  l2_err = fem_l2_norm(e_h);
  h10_err = fem_h10_norm(e_h);
  disp(['L2-error norm |u_h - I(u_exact)| = ',num2str(l2_err)]);
  disp(['H10-error norm |u_h - I(u_exact)| = ',num2str(h10_err)]);
  
  if ~no_plot
    disp('solution and errors plotted. Inspect workspace.')
    keyboard;
  end;
  
  res.l2_error = l2_err;
  res.h10_error = h10_err;
  
 case 4  % error convergence investigation
  
  params = [];
  params.dimrange = 1;
  pdeg = 2;
  qdeg = 3 * pdeg; % quadrature_degree
  params.pdeg = pdeg;
  params.qdeg = qdeg;
  params.no_plot = 1;
  
  numintervals = [2,4,8,16,32,64,128,256];
  l2_errs = zeros(length(numintervals),1);
  h10_errs = zeros(length(numintervals),1);
  for i = 1:length(numintervals)
    disp(['----------------------------------------------']);
    disp(['numintervals = ',num2str(numintervals(i))]);
    params.numintervals = numintervals(i);
    res = fem_poisson(3,params);
    l2_errs(i) = res.l2_error;
    h10_errs(i) = res.h10_error;
  end
     
 case 5 % detailed and rb simulation by high-level methods
  
  disp('detailed simulation:');
  model = elliptic_debug_model;
  model = elliptic_discrete_model(model);
  model_data = gen_model_data(model);
  model = model.set_mu(model,[0,0,0]);
  sim_data = detailed_simulation(model,model_data);
  plot_params.title = 'detailed_solution';
  figure, plot_sim_data(model,model_data,sim_data,plot_params);
 
  disp('offline computations (gen_detailed_data: basis):');
  detailed_data = gen_detailed_data(model,model_data);
  disp('offline computations II (gen_reduced_data: operator components):');
  reduced_data = gen_reduced_data(model,detailed_data);

  disp('online reduced simulation:');
  model = model.set_mu(model,[0,0,0]);
  rb_sim_data = rb_simulation(model,reduced_data);
  rb_sim_data = rb_reconstruction(model, detailed_data, rb_sim_data);
  plot_params.title = 'reduced solution';
  figure, plot_sim_data(model,model_data,rb_sim_data,plot_params);
  
  eh = sim_data.uh - rb_sim_data.uh;
  l2_err = fem_l2_norm(eh);  
  h10_err = fem_h10_norm(eh);  
  disp(['L2-error: ',num2str(l2_err)]);
  disp(['H10-error: ',num2str(h10_err)]);

 case 6 
  
  disp('demo_rb_gui:')  
  model = elliptic_debug_model;
  model = elliptic_discrete_model(model);
  model_data = gen_model_data(model);
  detailed_data = gen_detailed_data(model,model_data);
  plot_params.axis_tight = 1;
  demo_rb_gui(model,detailed_data,[],plot_params);
    
end;