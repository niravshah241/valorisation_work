function p = ldg_plot(dofs,grid,params)
%function p = ldg_plot(dofs,grid,params)
%
% function plotting a single scalar ldg function on triangular grid
% A patch plot is performed as default using speified number of
% subsamplings of the triangles. Each patch is plotted linearly 
% interpolated (hence only showing true values in subsampling nodes.)
% On each triangle-edge, subsamp_level points are inserted.
%
% p is the list of handles to the graphics primitives
%
% grid must provide the fields 
%      X,Y,VI,CX,CY,nelements,nvertices and nneigh;
%
% params must specify the ldg-type: params.dimrange, params.pdeg
%
% optional fields of params:
%   shrink_factor : if this flag is given, the elements are plotted shrinked
%   axis_equal    : if this flag is set, set axis to equal scale
%   no_lines      : if this flag is set, no lines are drawn
%   show_colorbar : if this flag is set, a colorbar is drawn (default 1)
%   colorbar_location : string specifying the position of the
%                 colorbar, e.g. 'South','EastOutside' (default), etc.
%   clim          : if this 2-vector is set, the colorbar is set to
%                   these values

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


% Bernard Haasdonk 2.2.2009

if nargin<3
  params = [];
end;

if isnumeric(dofs)
  df = ldgdiscfunc(dofs,params);
else
  df = dofs;
end;
df.grid = grid;
p=plot_discfunc(df,params);

