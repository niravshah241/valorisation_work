function einds_ext = index_ext(grid, einds, neigh_steps)
%function einds_ext = index_ext(grid, einds, neigh_steps)
%
% function computing the sorted list of element indicés, which are
% given by the neigh_steps number of neighbours of the elements with
% indices einds in grid.
% e.g. for neigh_steps = 0: einds_ext = einds 
% e.g. for neigh_steps = 1: einds_ext = einds plus all edge-neighbours
%      of elements 
% e.g. for neigh_steps = 2: einds_ext = einds plus all edge-neighbours
%      of elements plus all edge neighbours of these.

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

mask = zeros(1,grid.nelements);
mask(einds) = 1;
einds_ext = einds;
for n = 1:neigh_steps
  nbi  = grid.NBI(einds_ext,:);
  i    = find(nbi>0);
  mask(nbi(i)) = n+1;
  einds_ext = find(mask>0);
end;