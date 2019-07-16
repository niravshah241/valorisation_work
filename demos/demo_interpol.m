function demo_interpol
% function femo_interpol
%
% Script using femdiscfunc for interpolating functions by local and
% global evaluations.

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

disp('-------------------------------');
disp('     FE-interpolation demo     ');
disp('-------------------------------');

% poisson_model
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

%df = fem_interpol_local(model.source, df);
df = fem_interpol_global(model.source, df); 
figure,plot(df);
title('source function');
