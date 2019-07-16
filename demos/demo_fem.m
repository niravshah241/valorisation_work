function res = demo_fem(params)
% function res = demo_fem(params)
%
% script solving fem-problem, plotting solution and errors,
% return values: errors.

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


% B. Haasdonk, I. Maier 26.04.2011

% finite element assembly, solution and plot.
% different components assembled separately, as later
% affine-decomposition will require such modularization

% assemble discrete system, solve and plot results

disp('------------------------');
disp('        FEM demo        ');
disp('------------------------');

if nargin == 0
  params = [];
end;

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
model = elliptic_discrete_model(model);
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
  
res.l2_error = l2_err;
res.h10_error = h10_err;