function edges = lagrange_nodes_edges(pdeg);
%function edges = lagrange_nodes_edges(pdeg);
% 
% function returning a matrix indicating the containment of
% lagrange nodes in local edges.
% edges(lid,eid) 
% => 1 if lagrange node with index lid is lying on edge local index eid
% => 0 if lagrange node with index lid is not lying on edge local eid

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


% Bernard Haasdonk 12.1.2011

if pdeg <= 0
  error('pdeg not supported!');
end;

switch pdeg
 case 1
  edges = [1, 0, 1; 1, 1, 0; 0, 1, 1];
 case 2
  edges = [1, 0, 1; 1, 0, 0; 1, 1, 0; 0, 0, 1; 0, 1, 0; 0, 1, 1];
 case 3
  edges = [1, 0, 1; 1, 0, 0; 1, 0, 0; 1, 1, 0; 0, 0, 1; 0, 0, 0; ...
	   0, 1, 0; 0, 0, 1; 0, 1, 0; 0, 1, 1];
 otherwise
  lagrange_nodes = lagrange_nodes_lcoord(pdeg);
  edges = zeros(size(lagrange_nodes,1),3);
  i = find(abs(lagrange_nodes(:,2))<=eps); % i.e. lower edge 1
  edges(i,1)= 1;  
  i = find(abs(lagrange_nodes(:,1))<=eps); % i.e. left edge 3
  edges(i,3)= 1;  
  i = find(abs(1-lagrange_nodes(:,1)-lagrange_nodes(:,2))<=eps); % i.e. right edge 2
  edges(i,2)= 1;  
end;