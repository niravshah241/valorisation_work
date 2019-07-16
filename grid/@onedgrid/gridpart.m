function gridp = gridpart(grid,eind)
%function gridp = gridpart(grid,eind)
%
% function extracting a part of a triagrid, rectgrid or onedgrid 
% defined by the given 
% element indices in the vector eind. The neighbour information of
% the new resulting boundaries is set to -10

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


% Bernard Haasdonk 18.6.2010

%keyboard;
gridp = onedgrid(grid);
gridp.X = grid.X(eind);
gridp.nelements = length(eind);
gridp.global_eind = grid.global_eind(eind); % store original element indices

% generate elementid translation map: T \mapsto T_{local} 
new_el_id = zeros(1,grid.nelements);
new_el_id(eind) = 1:length(eind);

gridp.NBI = grid.NBI(eind,:);
i = find(gridp.NBI>0);
gridp.NBI(i) = new_el_id(gridp.NBI(i));
i = find(gridp.NBI == 0);
if ~isempty(i)
  gridp.NBI(i)= -10;
end;

