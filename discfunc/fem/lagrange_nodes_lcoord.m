function lcoords = lagrange_nodes_lcoord(pdeg);
%function lcoords = lagrange_nodes_lcoors(pdeg);
%
% function returning the local coordinates of the Lagrange nodes of
% the reference triangle. Result is a numlagrangepoints x 3 matrix

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
  lcoords = [0,0; 1,0; 0,1];
 case 2
  lcoords = [0,0; 0.5,0; 1,0; 0,0.5; 0.5,0.5; 0,1];
 case 3
  lcoords = [0,0; 1,0; 2,0; 3,0; 0,1; 1,1; 2,1; ...
	         0,2; 1,2; 0,3]/3;
 otherwise
%  lcoords = zeros((pdeg+1)*(pdeg+2)/2,3);
  c1 = (0:pdeg)'*ones(1,pdeg+1);
  c2 = c1';
  c3 = pdeg - c1 - c2;
  i = find(c3>=-eps);
  lcoords = [c1(i),c2(i)]/pdeg;
end;
