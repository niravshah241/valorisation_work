function df = fem_interpol_local(f, df, params)
%function df = fem_interpol_local(f, df, params)
%
% function interpolating the given analytical function with 
% local evaluation f(grid,elids,lcoord,params) into the allocated discrete
% function df
% by this, also discrete functions on the same grid with possibly
% other polynomial degree or other discrete function type can 
% be interpolated.

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


% B. Haasdonk 12.1.2011

% for all lagrange nodes evaluate f and set as dof

if nargin<3
  params = [];
end;

lagrange_nodes = lagrange_nodes_lcoord(df.pdeg); 

% possibly multiple evaluation in edge nodes, accept for now. 
for i = 1:size(lagrange_nodes,1);
  gids = df.global_dof_index(1:df.grid.nelements, i);
%  glob = local2global(df.grid, 1:df.grid.nelements, lagrange_nodes(i,:));
  fs = f(1:df.grid.nelements,lagrange_nodes(i,:),params);
  df.dofs(gids) = fs;
end;

