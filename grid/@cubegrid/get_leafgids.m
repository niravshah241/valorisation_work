function gids = get_leafgids(grid)
%function gids = get_leafgids(grid)
% return leaf ids of cubegrid, i.e. vector with global element indices of leaf
% elements
%
% Return values:
%  vector of global leaf indices
%
% Bernard Haasdonk 27.3.2007

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

  
gids = find(grid.isleaf);  

end
