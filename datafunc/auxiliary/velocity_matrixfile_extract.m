function [Vx, Vy, lambda] = velocity_matrixfile_extract(fullfn,X,Y)
%function [Vx, Vy, lambda] = velocity_matrixfile_extract(fullfn,X,Y)
%
% function extracting the velocities from the given file in the
% points X and Y. A search is performed in case of partial access. 
% So this function is expensive. 

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

  
%V = load(fullfn);
V = cache('load',fullfn);

% find integer indices such that V.X(j)==X(i) and V.Y(j) = Y(i) 
% if whole flux-matrix is requested

if ~isempty(X) % otherwise size 0/0 => 0/1 !!!
  X = X(:); Y = Y(:);
end;

if isequal(X,V.X(:)) && isequal(Y,V.Y(:))
  i = 1:length(X);
  j = 1:length(X);
elseif isempty(X)
  i = [];
  j = [];
else
  % otherwise search corresponding matrix entries
  %% the following gives memory problem in 'large point set' requests
  %	keyboard;
  %	VXX = repmat(V.X(:)',length(X),1);
  %	VYY = repmat(V.Y(:)',length(X),1);
  %	XX = repmat(X(:),1,length(V.X));
  %	YY = repmat(Y(:),1,length(V.X));
  %	[i,j] = find(VXX==XX & VYY==YY); 
  % should give one entry for each i!!
  % if ~isequal(sort(i),(1:length(X))')
  %  error('requested points not or multiple times in matrix!');
  % end;
  
  % so mixture: blockwise approach as above
  %	i = [];j = [];
  %	% compute optimal cut-size if required
  %	n1 = length(X(:));
  %	cutsize = ceil(5000000/length(V.X));
  %	ncuts = ceil(n1/cutsize)
  %	for c = 1:ncuts
  %	  ind1 = (1+(c-1)*cutsize):min(n1, c*cutsize);
  %	  ind1 = ind1';
  %	  VXX = repmat(V.X(:)',length(ind1),1);
  %	  VYY = repmat(V.Y(:)',length(ind1),1);
  %	  XX = repmat(X(ind1),1,length(V.X));
  %	  YY = repmat(Y(ind1),1,length(V.X));
  %	  [ipart,jpart] = find(VXX==XX & VYY==YY);
  %	  j = [j; jpart(:)];
  %	  i = [i; ipart(:)+(c-1)*cutsize];
  %	end;	
  
  [i,j] = find_corresp([X(:)';Y(:)'],[V.X(:)';V.Y(:)']);
  % remove occasional double find-results (inner edges)
  [si, sind] = sort(i);
  sj = j(sind);
  % forward-difference 
  diff = si(2:end)-si(1:end-1);
  k = find(diff~=0); % => si(k+1)~=si(k) 
  i = [si(1);si(k+1)];
  j = [sj(1);sj(k+1)];
  %should give one entry for each i!!
  if ~isequal(sort(i),(1:length(X))')
    error('requested points not or multiple times in matrix!');
  end;
end;
% i(k) is the linear index of the point k in the list X / Y 
% j(k) is the linear index of the point k in the list V.X / V.Y 
Vx = NaN*ones(size(X));
Vy = NaN*ones(size(X));
if ~isempty(i)
  Vx(i) = V.Vx(j);
  Vy(i) = V.Vy(j);
end;
lambda = V.lambda;

%| \docupdate 
