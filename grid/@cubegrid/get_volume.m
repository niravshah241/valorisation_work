function vols = get_volume(grid,gids)
% function vols = get_volume(grid,gids)
% method determining the volumes of a set of elements
%
% A simple product of edge-lengths is computed, i.e. the edges are assumed to
% be axis-parallel.
%
% Parameters:
%  gids: Global indices
%
% Return values:
%  vols: The volumes for the elements given by 'gids'.
%
% @note The computation for multiple gids is trivially performed by a loop.
% Vectorization must be done sometime!!

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
  
  vols = ones(size(gids));
  dim = grid.dimension;
  
  % determine edge indices, which are to be used:
  % 0D: eid(0) = [] of ne(0) = 0 edges 
  % 1D. eid = [1] of ne(1) = 1 edges
  % 2D. eid = [1, 3] of ne(2) = 4 edges
  % 3D. eid = [1, 3, 9] of ne(3) = 12 edges 
  % ND: eid(N) = [eid(N-1), 2* ne(N-1)+1] of ne(N) = 2*ne(N-1) + 2^(N-1) edges
  
  eid = [];
  ne = 0;
  for i = 1:dim
    eid = [eid, 2*ne+1];
    ne = 2* ne + 2^(i-1);
  end;
  
  for gid = 1:length(gids);
    [p0,p1] = get_edges(grid,gids(gid));
    diff = p0(:,eid)-p1(:,eid);
    len = sqrt(sum(diff.^2));
    vols(gid) = prod(len);     
  end;
end

