function indices = regular_submatrix_indices(M,epsilon)
%function indices = regular_submatrix:indices(M,epsilon)
%
% function determining the indices of a submatrix of the square
% matrix M, such that M(indices,indices) is a regular matrix with 
% determinant larger
% than epsilon. For this, simple matrix-growing is performed
% starting with the left upper element. If 
% addition of a row/column leads to too-small determinant, this is skipped.
%
% this routine can be applied on gram-matrices to detect linear
% independent vectors.

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


% Bernard Haasdonk 14.11.2007

indices = 1;

for j = 2:size(M,1);
  ind_new = [indices, j];
  B = M(ind_new,ind_new);
  if abs(det(B))>=epsilon
    % add j row/column
    indices = ind_new;
  end
end

% TO BE ADJUSTED TO NEW SYNTAX
%| \docupdate 
