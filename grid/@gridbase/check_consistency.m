function res = check_consistency(grid)
%function res = check_consistency(grid)
% function checking the consistency of a polygonal grid, i.e. checking, whether
% the edge centroids correspond (implicit neighbour-index test)
%
% @todo: perhaps later further extensions
%
% return values:
%  res: boolean value which is 'true' if the check succeeds, or 'false' if it
%  fails.
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


% Bernard Haasdonk 9.5.2007

  res = 1;

% check neighbour indices corresponding across edges
% must satisfy NBI(NBI(i,j), INB(i,j)) = i;
% works for general grids, but uses loops... 
%disp('checking neighbour-edge relations');
%disp('and neighbour-edge-centroid consistency:');
for i=1:grid.nelements
  for di=1:grid.nneigh
    if (grid.NBI(i,di)>0)
      if grid.NBI(grid.NBI(i,di), grid.INB(i,di))~=i
	disp('neighbour relations are inconsistent!!');
	res = 0;
      end;
      if (abs(grid.ECX(grid.NBI(i,di), grid.INB(i,di))-grid.ECX(i,di))>eps)
	disp('edge ECX are inconsistent!!');
	res = 0;
      end;
      if (abs(grid.ECY(grid.NBI(i,di), grid.INB(i,di))-grid.ECY(i,di))>eps)
	disp('edge ECY are inconsistent!!');
	res = 0;
      end;
    end;
  end,
end;

