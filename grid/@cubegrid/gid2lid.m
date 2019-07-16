function lids = gid2lid(grid,gids)
%function lids = gid2lid(grid,gids)
% function converting global element ids to leaf element ids
%
% Parameters:
%  gids: global indices of elements to which the local indices shall be
%  computed.
%
% Return values:
%  lids: leaf element ids or '-1' if a global id was not a leaf element

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


% Bernard Haasdonk 27.3.2007
  
gids=gids(:);  
leafgids = find(grid.isleaf);  
lid_array = -1*ones(1,grid.nelements);
lid_array(leafgids) = 1:length(leafgids);
lids = lid_array(gids);
lids = lids(:);

end
