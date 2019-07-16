function display(grid)
%function display(grid)
% cubegrid's display method.
%
% All scalar quantities are displayed. Internal pointlists, vectors and flags
% of elements are not displayed.  These can be accessed by the get() method.

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
    
disp(['dimension : ', num2str(grid.dimension)]);  
rangestr = '';
%for i=1:grid.dimension
%  rangestr = [rangestr,' [ ',  num2str(grid.range{i}),' ]'];
%end;
%disp(['ranges : ', rangestr]);
%disp(['number of intervals : ',num2str(grid.numintervals)]);
nlevels = max(grid.level);
disp(['number of refinement steps : ',num2str(grid.refine_steps)]);
disp(['number of levels : ',num2str(nlevels)]);
for l=0:nlevels
  disp(['number of elements on level ',num2str(l), ' : ', ...
	num2str(length(find(grid.level==l)))]);
end;
disp(['number of leaf-elements : ', ...
      num2str(length(find(grid.isleaf)))]);
disp(['number of elements : ', num2str(grid.nelements)]);
disp(['number of vertices : ', num2str(grid.nvertices)]);


%display(grid.vertex);
%display(grid.vertexindex);
%display(grid.level);
%display(grid.isleaf);

end

