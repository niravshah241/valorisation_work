function ret = get(grid,propertyname, varargin)
%function ret = get(grid,propertyname, varargin)
% get-method for cubegrids
%
% Parameters:
%  propertyname: A property name of the grid which can be one of the following
%       - 'level' : get vector of levels of all elements
%       - 'dimension' : get dimension of grid
%       - 'nvertices' : get number of vertices of grid
%       - 'nelements' : get number of elements of grid
%       - 'isleaf' : get vector of isleaf-indices of all elements
%       - 'vertexindex' returns a matrix '(i,j)' index of 'j'-th vertex of element 'i'
%                     a further optional argument with element indices can be passed,
%                     then only the vertices of these elements are returned,
%                     @code get(grid,'vertexindex',[2,3]) @endcode returns a
%                       '2 x nvertices_per_element' matrix with indices
%       - 'vertex' : returns a matrix of all vertex coordinates, i.e. 
%                ret(i,j) => j-th coordinate of i-th vertex point
%                optionally, a subset of vertex-indices can be determined as
%                further parameter, i.e.
%                @code get(grid,'vertex',[2 3 4]) @endcode produces a '3 x dim'
%                matrix with the coordinates of the points
%       .
%       or property can be equal to one of the following strings returning
%       quantities which are derived from the above grid properties
%       - 'nleafelements' : get number of leaf elements of grid
%       - 'leafelements'  : get vector of indices of elements that are leaf elements
%       - 'leafcogs' : returns a matrix of all centers of gravity of the leaf
%              elements. 'ret(i,j)' is the 'j'-th coordinate of the 'i'-th leaf
%              element
%  varargin: An optional list of arguments as described above.

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


% Bernard Haasdonk 1.3.2007

  if nargin <2
    help('cubegrid/get');
    return;
  end;
  
  switch lower(propertyname)
   case 'level'
    ret = grid.level;
   case 'dimension'
    ret = grid.dimension;
   case 'nvertices'
    ret = grid.nvertices;
   case 'nelements'
    ret = grid.nelements;
   case 'nleafelements'
    ret = sum(grid.isleaf);
   case 'isleaf'
    ret = grid.isleaf;
   case 'range'
    ret = grid.range;
   case 'numintervals'
    ret = grid.numintervals;
   case 'vertexindex'
    if nargin <= 2
      ids = 1:grid.nelements;
    else
      ids = varargin{1};
    end;
    ret = grid.vertexindex(ids,:);
   case 'vertex'
    if nargin <= 2
      ids = 1:grid.nvertices;
    else
      ids = varargin{1};
    end;
    ret = grid.vertex(ids,:);
   case 'leafelements'
    ret = find(grid.isleaf);
   case 'leafcogs'
    nleafelements = sum(grid.isleaf);
    ids = 1:nleafelements;    
    ret = zeros(nleafelements,grid.dimension);
    leafvind = grid.vertexindex(find(grid.isleaf),:); 
    leafv = leafvind(:);
    % fill all components of ret:
    for i = 1:size(ret,2)
      C = grid.vertex(leafv,i);
      CC = reshape(C,size(leafvind));
      ret(:,i) = mean(CC,2);
    end;   
   otherwise
    error('Propertyname not supported!');
  end;

end

