function ngrid = remove_duplicate_vertices(grid , epsilon )
%function ngrid = remove_duplicate_vertices(grid [, epsilon])
% method to be used, if a vertex list in the grid should be compressed
%
% This might be required after refinement of a grid, merging of grids, etc.
% Duplicate detection is based on l2-error deviation thresholded by 'epsilon'.
%
% Parameters:
%  epsilon: l2-error deviation threshold. (Default = '1e-6')
%
% Return values:
%  ngrid: the compressed grid

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


% Bernard Haasdonk 15.3.2007
  
% detect duplicates in current vertex list wrt. epsilon derivation
  
  if nargin < 2
    epsilon = 1e-6;
  end;
  
  [un_id, dupl_id, un_id_of_dupl] = ...
  detect_duplicate_rows(grid.vertex, epsilon);
  % i.e. vector with index dupl_id(1) is identical to vector un_id_of_dupl(1)
    
  ngrid = grid;  

  % only perform compression, if required
  if ~isempty(dupl_id)
    % generate compressed vertex list  
    ngrid.vertex = ngrid.vertex(un_id,:);

    % construct vertex-index-translation map 
    % such that new_id(dupl_id(1)) = un_id_of_dupl(1)
    % and new_id(un_id(1)) = un_id(1);
    new_id = -ones(max(max(grid.vertexindex)),1);
    new_id(un_id) = 1:length(un_id);
    new_id(dupl_id) = new_id(un_id_of_dupl);

    % perform vertex-index-translation
    vertexindex_vector = new_id(grid.vertexindex(:));      
    ngrid.vertexindex = reshape(vertexindex_vector, ...
				 size(grid.vertexindex));
    
    if ~isempty(find(vertexindex_vector <= 0,1))
      error('translation map erroneous');
    end;

    ngrid.nvertices = size(ngrid.vertex,1);
    
  end;

end

