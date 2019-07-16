function gids = lid2gid(grid,lids)
%function gids = lid2gid(grid,lids)
% function converting leaf-element ids to global element ids
%
% Parameters:
%  lids: local indices for which the global ones shall be computed
%
% Return values:
%  gids: global indices of leaf indices

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

lids = lids(:);  
leafgids = find(grid.isleaf);  
gids = leafgids(lids);

end

