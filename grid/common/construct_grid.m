function grid = construct_grid(params)
%function grid = construct_grid(params)
%
% function giving a simple common constructor syntax to different
% grid types. See the constructors of the grid types with a single
% params option for details on specifications
%
% required fields of params:
%    gridtype : 'rectgrid' or 'triagrid'

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


% Bernard Haasdonk 18.5.2007

if isfield(params,'gridtype')
  switch params.gridtype
   case 'rectgrid'
    grid = rectgrid(params);
   case 'onedgrid'
    grid = onedgrid(params);
   case 'triagrid'
    grid = triagrid(params);
    grid = set_boundary_types(grid, params);
   case 'none'
    % nothing to do
    grid = [];
   otherwise
    error('gridtype not known')
  end;
else
  grid = [];
end;
