function p = fem_plot_dofmap(df,params)
%function p = fem_plot_dofmap(df,params)
%
% function plotting the numbers of the dofs for optical consistency
% check across edges

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

p = plot(df.grid);
hold on
lagrange_nodes = lagrange_nodes_lcoord(df.pdeg);
for li = 1:size(lagrange_nodes,1);
  clocal = repmat([1/3,1/3],size(lagrange_nodes,1),1);
  % shrink location a bit towards center for readability;
  local = (lagrange_nodes - clocal)*0.8 + clocal;
  glob = local2global(df.grid,1:df.grid.nelements,local(li,:));
  nums = num2str(df.global_dof_index(:,li));
  p2 = text(glob(:,1),glob(:,2),nums,'Fontsize',5);
end;
p = [p;p2];