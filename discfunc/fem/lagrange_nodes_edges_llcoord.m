function llcoord = lagrange_nodes_edges_llcoord(pdeg)
% function llcoord = lagrange_nodes_edges_llcoord(pdeg)
%
% function returning edge-local coordinates of the
% lagrange-nodes of the reference triange. llcoord is a
% nlagrange_nodes x 3 matrix.
% llcoord(i,j) = real number on the unit-interval, if the i-th
%                lagrange-node lies on the j-th edge
% llcoord(i,j) = -1, else

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


% Immanuel Maier, 11.03.2011

switch pdeg
 case 1
  llcoord = [1 -1 0;
	     0 1 -1;
	     -1 0 1];
 case 2
  llcoord = [1 -1 0;
	     0.5 -1 -1;
	     0 1 -1;
	     -1 -1 0.5;
	     -1 0.5 -1;
	     -1 0 1];
 case 3
  llcoord = [3 -3 0;
	     2 -3 -3;
	     1 -3 -3;
	     0 3 -3;
	     -3 -3 2;
	     -3 -3 -3;
	     -3 2 -3;
	     -3 -3 1;
	     -3 1 -3;
	     -3 0 3]/3;
 otherwise
  error('pdeg < 1 or > 3 not supported in lagrange_nodes_edges_llcoord!');
end;