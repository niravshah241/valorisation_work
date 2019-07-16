function res = isequal(grid,grid2)
%function res = isequal(grid,grid2)
% function comparing two grids field by field
%
% parameters:
%  grid2: the grid we want to compare with
%
% Return values:
%  res: boolean indicating whether the two grids are equal

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


% Bernard Haasdonk 15.2.2011

fields = {'nelements','nvertices','nneigh','A','Ainv','VI','X','Y','CX',...
	  'CY','NBI','INB','EL','DC','NX','NY','ECX','ECY','SX','SY','ESX',...
	  'DS','hmin','alpha','JIT'};

res = 1;
for i = 1:length(fields)
  res = res && isequal(getfield(grid,fields{i}),getfield(grid2,fields{i}));
end

