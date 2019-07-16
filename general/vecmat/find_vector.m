function ind = find_vector(vec, mat, epsilon)
%function ind = find_vector(vec, mat[, epsilon])
%
% find a vector in a matrix. If vec is a row vector, rows of mat are
% searched, if vec is a column, columns of mat are searched.
% equality is decided by l2-error being less than epsilon. epsilon can be
% ommited, is default 1e-6

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

      
% Bernard Haasdonk 29.3.2007
  
% if row vector, transpose the problem
  
  if nargin < 3
    epsilon = 1e-6;
  end;
  
  if size(vec,1) == 1
    vec = vec';
    mat = mat';    
  end;
  
  vvec = repmat(vec,1,size(mat,2));
  diff = sqrt(sum((vvec-mat).^2));
  ind = find(diff<epsilon);
    
  
% TO BE ADJUSTED TO NEW SYNTAX
%| \docupdate 
