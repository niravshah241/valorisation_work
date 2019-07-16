function gridp = gridpart(grid,eind)
%function gridp = gridpart(grid,eind)
% function extracting a part of a ::triagrid or ::rectgrid defined by the given
% element indices in the vector 'eind'.
%
% @note The neighbour information of the new resulting boundaries is set to
% '-10'
%
% @note The properties gridbase.hmin, gridbase.alpha and the
% distance-information in the new boundary elements are simply copied. I.e.
% these fields do not completely meet the definition in the constructor. They
% might be chosen slightly different, such that the 'gridp' would be really
% identical to the result generated by the constructor on the subset of points.
%
% Parameters:
%  eind: vector of cell indices which shall be extracted from the grid.
%
% Return values:
%  gridp: the partial grid of type ::gridbase with extracted cells 'eind'.

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


% Bernard Haasdonk 15.5.2007

gridp = copy(grid);
gridp.nelements = length(eind);

% generate elementid translation map: T \mapsto T_{local} 
new_el_id = zeros(1,grid.nelements);
new_el_id(eind) = 1:length(eind);

% generate vertex translation map: V \mapsto V_{local}
new_vertex_id = zeros(1,grid.nvertices);
mask = zeros(1,grid.nvertices);
mask(grid.VI(eind,:))= 1;
vind = find(mask);
new_vertex_id(vind) = 1:length(vind);

% set number of vertices, area vector, and inverse area vector
gridp.nvertices = max(new_vertex_id);
gridp.A = gridp.A(eind);
gridp.Ainv = gridp.Ainv(eind);

% set vertex coordinates
gridp.X = grid.X(vind);
gridp.Y = grid.Y(vind);

% set VI map
gridp.VI = grid.VI(eind,:); %: VI(i,j) is the global vertex index of j-th
gridp.VI = new_vertex_id(gridp.VI);

% set coordinates of elements' barycenters
gridp.CX = gridp.CX(eind);
gridp.CY = gridp.CY(eind);

% ATTENTION! Set neighbour indices
gridp.NBI = grid.NBI(eind,:);
i = find(gridp.NBI>0);
gridp.NBI(i) = new_el_id(gridp.NBI(i));
i = find(gridp.NBI == 0);
if ~isempty(i)
  gridp.NBI(i)= -10;
end;

gridp.INB = grid.INB(eind,:);

gridp.EL = grid.EL(eind,:);
gridp.DC = grid.DC(eind,:);
gridp.NX = grid.NX(eind,:);
gridp.NY = grid.NY(eind,:);
gridp.ECX = grid.ECX(eind,:);
gridp.ECY = grid.ECY(eind,:);
gridp.SX = grid.SX(eind,:);
gridp.SY = grid.SY(eind,:);
gridp.ESX = grid.ESX(eind,:);
gridp.ESY = grid.ESY(eind,:);
gridp.DS = grid.DS(eind,:);

