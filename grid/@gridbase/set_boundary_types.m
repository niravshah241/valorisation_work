function grid = set_boundary_types(grid,params)
%function grid = set_boundary_types(grid,params)
% function setting the boundary types of a polygonal grid.
%
% For ::rectgrid's this is already done in the constructor. For ::triagrid this
% can be performed explicitly here. Existing boundary settings are 'painted' by
% rectangles The edges with midpoints within such a rectangle are marked
% accordingly.
%
% required fields of params:
%    bnd_rect_corner1: coordinates of lower corner of to be marked
%              boundaries. This can also be cell array of coordinates for
%              drawing different rectangles.
%    bnd_rect_corner2: coordinates of upper corner of to be marked
%              boundaries. This can also be cell array of coordinates for
%              drawing different rectangles.
%    bnd_rect_index: integer index to be set on the edges in the above  defined
%            rectangle. Should not be positive integer in the range of the
%            number of elements. Use negative indices for certain later
%            discrimination.
%              - '-1' is treated as Dirichlet in the numerics
%              - '-2' is treated as Neumann in the numerics
%              .
%            This can also be a vector of the same length as the cell arrays
%            defined for 'bnd_rect_corner1' and 'bnd_rect_corner2' in case of
%            different boundary types on the domain.
%

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


% Bernard Haasdonk 14.5.2007

if isfield(params,'bnd_rect_index')
  %    keyboard;
  bnd_ind = find(grid.NBI<=0);  
  SX = grid.ECX(bnd_ind);
  SY = grid.ECY(bnd_ind);
  if (max(params.bnd_rect_index)>0)
    error('boundary indices must be negative!');
  end;
  if size(params.bnd_rect_corner1,1) == 1
    params.bnd_rect_corner1 = params.bnd_rect_corner1';
  end;
  if size(params.bnd_rect_corner2,1) == 1
    params.bnd_rect_corner2 = params.bnd_rect_corner2';
  end;
  for i = 1:length(params.bnd_rect_index)
    indx = (SX > params.bnd_rect_corner1(1,i)) & ...
	   (SX < params.bnd_rect_corner2(1,i)) & ... 
	   (SY > params.bnd_rect_corner1(2,i)) & ...
	   (SY < params.bnd_rect_corner2(2,i));
    grid = set_nbi(grid,bnd_ind(indx), ...
			params.bnd_rect_index(i));
  end;
end;

