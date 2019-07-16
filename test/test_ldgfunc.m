function OK = test_ldgfunc
% function testing ldgfunc abilities
%
% currently testing:
%   l2projection of vectorial function
%   extraction of scalar component
%   local evaluation of vectorial function
%   ldg_l2error

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


% Bernard Haasdonk 18.1.2010

OK = 1;

% initialize grid
grid = triagrid();

params.nelements = grid.nelements; 
params.pdeg = 1;
params.dimrange = 2; % vectorial function
params.ndofs = ldg_ndofs(params);
params.ndofs_per_element = ldg_ndofs_per_element(params);

% function constant [1.0, 1.0]
f = @(einds,loc,grid,params) ones(length(einds),1)*[1.0, 1.0];

%f = @(einds,loc,grid,params) f_global(...
%    local2global(grid,einds,loc,params),params); 

qdeg = 4;
dofs = ldg_l2project(f,qdeg,grid,params);
[dofs1, params1 ] = ldg_scalar_component(dofs,1,params);
%disp('local evaluation should be constant [1 1]');

res = ldg_evaluate(dofs,1:10,[0,0],grid,params);
if max(max(abs(res-ones(10,params.dimrange))))>1e-6
  disp('test_ldgfunc: local evaluation result of projected dofs not [1,1]!!');
  OK = 0;
end;

res = ldg_evaluate(dofs1,1:10,[0,0],grid,params1);
if max(max(abs(res-ones(10,1))))>1e-6
  disp('test_ldgfunc: local evaluation of scalar component not 1!!');
  OK = 0;
end;

df1 = ldgdiscfunc(params);
df1.dofs = dofs;
df1.grid = grid; %for now
res = df1(1:10,[0,0],params);
if max(max(abs(res-ones(10,params.dimrange))))>1e-6
  disp('test_ldgfunc: local evaluation result of @ldgdiscfunc not 1!!');
  OK = 0;
end;

res = ldg_l2error(dofs,f,qdeg,grid,params);
if abs(res)>1e-6
  disp('test_ldgfunc: ldg_l2error not 0!!');
  OK = 0;
end;

%disp('plot should be constant 1');
%ldg_plot(dofs1, grid, params1);
%keyboard;%| \docupdate 
