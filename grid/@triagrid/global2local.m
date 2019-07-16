function loc = global2local(grid, elementid, glob)
%function loc = global2local(grid, elementid, glob)
%
% function getting a triagrid, an element-ID and a vector of points and
% giving a vector of transformed points. The triangle given by elementid is
% used for the creation of an affine map to the standard tringle, then this
% transformation is used for all the points in glob
%
% input:
% grid:       triagrid
% elementid:  scalar nr of the triangle in triagrid, which defines the affine map
% glob:       n-by-2-vector of points to be transformed
%
% output:
% loc:        n-by-2-vector of the transformed points
%
% Oliver Zeeb, 02.02.11

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


vertex_id = grid.VI(elementid,:);
tria_pts_x = grid.X(vertex_id);
tria_pts_y = grid.Y(vertex_id);

[C,G] = aff_trafo_glob2loc(tria_pts_x, tria_pts_y);

loc = (repmat(C,1,size(glob,1)) + G*glob')';
