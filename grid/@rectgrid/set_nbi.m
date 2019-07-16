function grid = set_nbi(grid,nbind,values)
%function grid = set_nbi(grid,nbind,values)
% function setting some neighbour indices of a grid to specified values.
%
% Parameters:
%  nbind:  neighbour indices
%  values: new neighbour indices. This can be a single scalar or a vector of
%          the same length as 'nbind'.

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

grid.NBI(nbind) = values;

