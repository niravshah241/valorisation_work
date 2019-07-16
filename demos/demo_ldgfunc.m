function demo_ldgfunc
%function demo_ldgfunc
%
% demo for showing ldgfunc capabilities
% idea is mainly to keep dofs and size information separate in
% order to easily work with vectors and matrices of dofs for
% efficiency.
%
% A slower alternative can be beneficial: by using the \\@ldg function
% class, a dof-vector can be stored in the object. And hence an
% evaluation be performed by the subsref method as for analytical functions. 
% So ldg-objects can be used identically as any
% analytical function in integration, matrix assembly, etc.
% As classes/methods are slower than
% structs, this should only be used if necessary. 
% So most ldg operations will accept separate dof vectors/matrices
% and size information.

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


% Bernard Haasdonk 28.1.2009

% quadrature degree for integration:
qdeg = 4;
% initialize grid
grid = triagrid();

disp('');
disp('initialization of zero ldg function');

% initialize basic data
pdeg = 1;
dimrange = 2;
params = ldg_params(grid.nelements,pdeg,dimrange);
%params.nelements = grid.nelements; 
%params.pdeg = 1;
%params.dimrange = 2; % vectorial function
%params.ndofs = ldg_ndofs(params);
%params.ndofs_per_element = ldg_ndofs_per_element(params);
dofs = ldg_zero(params); % simple zero vector of dofs

% these variables are sufficient to perform most effective operations
% for ldg functions, i.e. data and parameters are kept disjoint.
%
% A slower alternative can be beneficial: by using the @ldg function
% class, a dof-vector can be stored in the object. And hence an
% evaluation be performed by the subsref method as for analytical functions. 
% So ldg-objects can be used identically as any
% analytical function in integration, matrix assembly, etc.
% As classes/methods are slower than
% structs, this should only be used if necessary:

df = ldgdiscfunc(dofs,params); % setting of data: persistent!
display(df);

% most of the following does not access df, but work with dofs and params

disp('----------------------------------------------------------');
disp(['dof-access and plot of scalar component with subsampling']);

% manipulate/read dofs and visualize
dofs = 1:params.ndofs;
[dofs_scalar1,sparams] = ldg_scalar_component(dofs,1,params);
for ssl = 0:2
  sparams.plot_subsampling_level = ssl;
  figure
  ldg_plot(dofs_scalar1,grid,sparams);
  title(['dofnumbers, subsampling level ',num2str(ssl)]);
end;
disp('press key to continue');
pause();

disp('');
disp('----------------------------------------------------------');
disp(['local evaluation of ldgfunc and analytical functions in point ',...
      ' lcoord=[1/3,1/3] in first 10 elements:']);
% define analytical function inline, e.g. constant 1
f1 = @(einds,loc,grid,params) ones(length(einds),1);
% define analytical function by pointer to matlab function 
f2 = @(einds,loc,grid,params) f_local(einds,loc,grid,params);
% define analytical function by pointer to matlab function with 
% internal conversion of local to global coordinates:
f3 = @(einds,loc,grid,params) f_global(...
    local2global(grid,einds,loc,params),params); 
% evaluate all functions
disp('f1 (scalar function constant 1):');
f1(1:10, [1/3,1/3], grid, params)
disp('f2: (function f_local defined in local coordinates)');
f2(1:10, [1/3,1/3], grid, params)
disp('f3: (function f_global defined in global coordinates)');
f3(1:10, [1/3,1/3], grid, params)
% note here a ldg function is used identically as the analytical functions!
disp('df: (@ldgdiscfunc dimrange=2, pdeg=1, zero-function)');
df.grid = grid;
df(1:10, [1/3,1/3], params) % identical call as above!
% this is equivalent to
% ldg_evaluate(dofs,1:10, [1/3,1/3], grid, params);

disp('press key to continue');
pause();

disp('');
disp('----------------------------------------------------------');
disp('test of orthogonality of element basis functions: gram matrix');

% test orthogonality of basis functions:
disp('gram matrix of basis on reference element should be identity:')
K = ldg_local_mass_matrix(qdeg,params)
%f = @(lcoord) gram_matrix(ldg_evaluate_basis(lcoord,params));
%triaquadrature(qdeg,f)

disp('press key to continue');
pause();

disp('');
disp('----------------------------------------------------------');
disp('l2 projection of analytic function: hor/vert sinus waves');

% l2 projection of an analytical function and plot
f = @(einds,loc,grid,params) f_global(...
    local2global(grid,einds,loc,params),params); 
dofs = ldg_l2project(f,qdeg,grid,params);
%dofs = l2project(f,qdeg,grid,params);

[dofs_scalar1,sparams] = ldg_scalar_component(dofs,1,params);
dofs_scalar2 = ldg_scalar_component(dofs,2,params);
sparams.plot_subsampling_level = 2;
figure,ldg_plot(dofs_scalar1,grid,sparams);title('component1');
figure,ldg_plot(dofs_scalar2,grid,sparams);title('component2');

disp('press key to continue');
pause();

disp('');
disp('----------------------------------------------------------');
disp('l2-error of zero ldg function with unity (should be 1)');

qdeg = 5;
pdeg = 3;
dimrange = 1;
params = ldg_params(grid.nelements,pdeg,dimrange)

% initialize zero function:
dofs0 = ldg_zero(params);
% initialize const-one function analytically 
f1 = @(einds,loc,grid,params) ones(length(einds),1);
% and project to ldgfunc
dofs1 = ldg_l2project(f1,qdeg,grid,params);

df1 = ldgdiscfunc(dofs1,params);
% compute errors with analytical and discrete function: 
% identical syntax by using ldg class!!!
params.evaluate = @ldg_evaluate;

err_df0_f1 = ldg_l2error(dofs0,f1,qdeg,grid,params);
disp(['l2 error (df0,f1) =',num2str(err_df0_f1)]);
err_df0_df1 = ldg_l2error(dofs0,df1,qdeg,grid,params);
disp(['l2 error (df0,df1) =',num2str(err_df0_df1)]);

disp('press key to continue');
pause();

disp('');
disp('----------------------------------------------------------');
disp(['Table of l2-errors of projection with increasing p, should' ...
      ' converge to 0']);
% l2-difference between ldgdiscfunc and function with 
% entity-local-coord evaluation possibility.
params.dimrange = 2;
params.nelements = grid.nelements;
qdeg = 10;
for p = 1:4
  params.pdeg = p;
  params.ndofs = ldg_ndofs(params);
  params.ndofs_per_element = ldg_ndofs_per_element(params);
  dofs = ldg_l2project(f,qdeg,grid,params);
%  dofs = l2project(f,qdeg,grid,params);
  l2err = ldg_l2error(dofs,f,qdeg,grid,params);
  disp(['pol-deg = ',num2str(p),', l2error(df,f) = ',num2str(l2err)]);
end;

disp('press key to continue');
pause();

return;

% l2-difference between two ldgdiscfuncs works identical
disp(['l2error(df,df) =',...
      num2str(l2error(df,df,qdeg,grid,params))]);

% arithmetics of ldg functions: sum, diff, scal-mtimes,
% in particular scalar product of vectorial localdgs



% right hand side assembly of LGS



% edge quadratures



% gradient operator



% div0 projection operator



% second operator



% concatenate for complete space operator




% implicit or explicit solution of time step




% stiffness-matrix assembly into sparse matrix




% solve FEM problem

% stop routine for inspection of workspace variables
keyboard;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function F = f_local(einds,loc,grid,params)
%function F = f_local(einds,loc,grid,params)
%
% function to be evaluated in local coordinates, i.e. 
% glob = [X, Y] with X and Y column vectors of global coordinates, 
% producing a corresponding
% result-matrix in F (each row one result)
% such functions can also be used with 3-dim coordinates without
% any change!
dimrange = 1;
F = zeros(size(einds,2),dimrange);
F(:,1) = einds;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function F = f_global(glob,params)
%function F = f_global(glob,params)
% function to be evaluated in gobal coordinates, i.e. 
% glob = [X, Y] with X and Y column vectors of global coordinates, 
% producing a corresponding
% result-matrix in F (each row one result)
% such functions can also be used with 3-dim coordinates without
% any change!
dimrange = 2;
F = zeros(size(glob,1),dimrange);
F(:,1) = sin(pi*glob(:,1));
F(:,2) = cos(pi*glob(:,2));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%| \docupdate 
