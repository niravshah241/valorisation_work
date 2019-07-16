function res = fem_ndofs(params,grid)
%function res = fem_ndofs_per_element(params,grid)
%
% function computing number of dofs based on 
% params.dimrange, params.grid and params.pdeg.
%
% Parameters:
%   grid:   object of type triagrid

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


% Bernard Haasdonk 2.9.2009


% num dofs on nodes:
ndof_nodes = grid.nvertices;

% num dofs on edges:
ndof_edges = (grid.nedges_interior + grid.nedges_boundary) * (params.pdeg-1);

% num dofs in interior:
ndof_elements = grid.nelements * (params.pdeg-1)* (params.pdeg-2)*0.5;

res = grid.nelements*params.dimrange*(params.pdeg+1)*(params.pdeg+2)/2;

res = (ndof_nodes + ndof_edges + ndof_elements) * params.dimrange;
