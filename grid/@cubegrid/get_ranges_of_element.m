function range = get_ranges_of_element(grid, element)
% return the ranges corresponding to the element of the
% grid, i.e. maximum vertex range of the part_num element
% of the grid.
%
% Markus Dihlmann 16.02.2010
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

if(element>get(grid,'nelements'))
    error('partition number bigger than actual number of leaf elements');
end


%get the vertices of this element
vertices = get(grid,'vertex',get(grid,'vertexindex',element));

dim=get(grid,'dimension');

range={};
for i=1:dim
    range{i}=[min(vertices(:,i)),max(vertices(:,i))];
end

end
