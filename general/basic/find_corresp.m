function [i,j]  = find_corresp(V1,V2)
%function [i,j]  = find_corresp(V1,V2)
%
% find identical columns in V1 and V2, i.e.
% V1(:,i(1)) = V2(:,j(1)) and all other indices in i,j
% double occurences can happen

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


% Bernard Haasdonk 22.5.2007

i = []; j = [];
n1 = size(V1,2);
n2 = size(V2,2);

if n1>0 && n2> 0
  cutsize = ceil(30000/n2);
  ncuts = ceil(n1/cutsize);
  for c = 1:ncuts
    ind1 = (1+(c-1)*cutsize):min(n1, c*cutsize);
    ind1 = ind1';
    mask = ones(length(ind1),n2);
    for dim = 1:size(V1,1);
      VV1 = repmat(V1(dim,ind1)',1,n2);
      VV2 = repmat(V2(dim,:),length(ind1),1);
      mask = mask & (VV1 == VV2);
    end;
    [ipart,jpart] = find(mask);
    j = [j; jpart(:)];
    i = [i; ipart(:)+(c-1)*cutsize];
  end;	
end;

%% for each column vector in the first matrix, the nearest neighbour
%% column in the second matrix is determined. such that
%% V(1) and V2(j(1)) are closest among all other V2 columns with
%% respect to l2-norm. The minimum distances are returned in dist
%%
%% The routine is required for large vector sets, such that the
%% simply tensor-product approach by setting up complete distance
%% matrices is avoided.
%%
%% Therefore, currently a stupid search is performed over all
%% columns of V1...
%
%% Bernard Haasdonk 22.5.2007
%
%n1 = size(V1,2);
%n2 = size(V2,2);
%j = zeros(n1,1);
%dmins = zeros(n1,1);
%for i=1:n1
%  dist = sum((repmat(V1(:,i),1,n2)-V2).^2); 
%  [dmin, jmin] = min(dist);
%  dmins(i) = sqrt(dmin(1));
%  j(i) = jmin(1);
%end;

% TO BE ADJUSTED TO NEW SYNTAX
%| \docupdate 
