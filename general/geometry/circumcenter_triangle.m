function c = circumcenter_triangle(q, p1, p2)
%function c = circumcenter_triangle(q, p1, p2)
%
% function computing the circumcenter of the triangle q,p1,p2.
%
% If q,p1,p2 is a matrix with columnwise points, then a matrix of 
% rownwise circumcenters is generated.

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
  error('Only 2d points acceptable in circumcenter computation');
end;


% let n1 be the normal to line q-p1
% let n2 be the normal to line q-p2
% then obviously c = (q+p1)/2 + lambda1 n1 = (q+p2) + lambda2 n2 
%
% let N = [n1 -n2]  
%
% => lambda = [lambda1; lambda 2] = N^-1 (p2-p1)/2 
%       where N^-1 = 1/Det(N) * [-n2y    n2x;  -n1y   n1x];
%
% so only one component lambda1 is required

n1 = [- (p1(:,2)-q(:,2)) , p1(:,1)-q(:,1)];
n2 = [- (p2(:,2)-q(:,2)) , p2(:,1)-q(:,1)];
detN = -n1(:,1).* n2(:,2) + n1(:,2).*n2(:,1);

lambda1 = 0.5 * detN.^(-1) .* (-n2(:,2).* (p2(:,1)-p1(:,1)) + ...
			      n2(:,1).* (p2(:,2)-p1(:,2)) );

c = 0.5 * (p1+q);
c(:,1) = c(:,1) + lambda1.* n1(:,1);
c(:,2) = c(:,2) + lambda1.* n1(:,2);

%| \docupdate 
