function leaf_element = coord2leaf_element(grid,coord)
%function coord2leaf_element(grid,coord)
%
%function giving the leaf_element of a grid where the values of the
%coord-vector can be found (coord- can for example be a set of mu-values in
%a parameter grid. This function utilizes a `O(log(n))`-search algorithm.
%
% coord :      vector of coord-values (must have same dimension as grid)
%
%
% 15.02.2010 Markus Dihlmann
%

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


if(get(grid,'dimension')~=length(coord))
    leaf_element=-1;
    error('length of coord-vector doesnt match the dimension of the grid')
end

dim = length(coord);
continue_search = 1;
element = 1;

% find maximal coordinates of vertexes to identify outer-bound coordinates
% of the whole cube.
V=get(grid,'vertex');
max=V(1,:);
for i=1:dim
  for j=1:length(V(:,1))
    if(max(i)<V(j,i))
      max(i)=V(j,i);
    end
  end
end

is_leaf = get(grid,'isleaf');

while continue_search
  leaf_element=element;
  ranges=get_ranges_of_element(grid,leaf_element);
  sum=0;
  for i=1:dim
    if coord(i)<max(i) %if coordinates are not on outer bound of grid
      if((coord(i)<ranges{i}(2))&&(coord(i)>=ranges{i}(1)))
        sum = sum+1;
      end
    else
      if((coord(i)<=ranges{i}(2))&&(coord(i)>=ranges{i}(1)))
        sum=sum+1;
      end
    end
  end
  if(sum == dim)
    if is_leaf(element)
      continue_search=0;
    else
      element = grid.firstchild(element);
    end
  else
    element=element+1;
  end

  if(element>get(grid,'nelements'))
    continue_search=0;
    leaf_element=-1;
    error('coord- range out of global parameter ranges');
  end
end

end


