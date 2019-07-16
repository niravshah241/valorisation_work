function a = area_triangle(q, p1, p2)
%function a = area_triangle(q, p1, p2)
%
% function computing the signed (!) area of the triangle q,p1,p2.
% is positive, if q,p1,p2 is counterclockwise, is negative, if
% q,p1,p2 is clockwise ordered.
%
% If q,p1,p2 is a matrix with rowwise points, then a vector of 
% areas is generated.

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
  error('Only 2d points acceptable in area computation');
end;

%a = zeros(size(p1,2),1);

% a = 0.5 * | Det (p1-q  ,  p2-q) |
a11 = p1(:,1)-q(:,1);
a21 = p1(:,2)-q(:,2);

a12 = p2(:,1)-q(:,1);
a22 = p2(:,2)-q(:,2);

a = 0.5 * (a11.*a22 - a21.*a12);
 


%| \docupdate 
