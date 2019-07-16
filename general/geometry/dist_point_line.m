function d = dist_point_line(q, p1, p2)
%function d = dist_point_line(q, p1, p2)
%
% function computing the distance of the 2d-point q from the line
% through the points p1, p2.
% If q,p1,p2 is a matrix with columnwise points, then a vector of 
% distances is generated.

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


% Bernard Haasdonk 10.5.2007


if size(p1,1)== 2
  p1 = p1';
end;

if size(p2,1)== 2
  p2 = p2';
end;

if size(q,1)== 2
  q = q';
end;

if (size(p1,2)~=2) || (size(p2,2)~=2) || (size(q,2)~=2)
  error('Only 2d points acceptable in distance computation');
end;

% d = area(triangle(q,p1,p2)*2 / dist(p1,p2))

% d = zeros(size(p1,2),1);
A = abs(area_triangle(q,p1,p2));
dp1p2 = p1 - p2;
d = A * 2 .* sqrt(dp1p2(:,1).^2 + dp1p2(:,2).^2).^(-1); 

%| \docupdate 
