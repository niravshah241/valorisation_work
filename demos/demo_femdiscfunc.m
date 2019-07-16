function demo_femdiscfunc
% function demo_femdiscfunc
%
% Script demonstrating some basic functionality of lagrange finite
% element functions.

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

disp('---------------------------------');
disp('      lagrange FE-functions      ');
disp('---------------------------------');

% poisson_model
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