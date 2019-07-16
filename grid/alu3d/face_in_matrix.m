function ind = face_in_matrix(face, matrix)
% function ind = face_in_matrix(face, matrix)
%
% function searching for the 4x1 vector face in the columns of the
% 4xn matrix by all cyclical permutations and reflection
% if found, the column-indices are returned in ind, otherwise it is empty

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

  
% Bernard Haasdonk 13.3.2006

matmat = repmat (matrix, 1,8);  % => 
ff = repmat(face(:),1,size(matrix,2));
faceface = [ff, ff([2 3 4 1],:), ff([3 4 1 2],:), ...
	     ff([4 1 2 3],:), ff([4 3 2 1],:), ff([3 2 1 4],:), ...
	     ff([2 1 4 3],:), ff([1 4 3 2],:)];
s = sum(abs(matmat-faceface));
i = find(s==0);
ind = [];
if ~isempty(i) 
  ind = mod(i-1,size(matrix,2))+1;
end;


% TO BE ADJUSTED TO NEW SYNTAX
%| \docupdate 
