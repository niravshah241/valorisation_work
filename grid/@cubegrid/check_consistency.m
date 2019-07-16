function res = check_consistency(grid)
%function res = check_consistency(grid)
%
% function checking the consistency of a cubegrid, i.e. checking, whether
% - all elements have axis parallel edges
% - the volume of the leaf-elements corresponds to the level-0 elements
%
% Return values:
%  res: This is 1 if OK, 0 if error
%
% @note The check is slow, as all edges for all elements are determined in
% loops and multiply

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


% perhaps later further extensions

% Bernard Haasdonk 27.3.2007

% for all elements check that all edges are axis-parallel  

  res = 1;
  
  for gid = 1:grid.nelements;
    [p0,p1] = get_edges(grid,gid);
    diff = p0-p1;
    % only one entry may relevantly differ from zero
    nonzeros = sum(abs(diff)>1e-10);
    if (min(nonzeros)<1) || (max(nonzeros)>1)
      disp('edge not axis parallel, please check!');
      res = 0;
    end;
  end;
  
  % check the volume of the leaf-elements corresponds to the level-0 elements
  
  macrogids = find(grid.level==0);
  leafgids = get_leafgids(grid);
  
  macrovol = get_volume(grid,macrogids);
  leafvol = get_volume(grid,leafgids);
  
  ratio = abs(sum(macrovol)-sum(leafvol))/sum(macrovol);
  if (ratio>1e-5)
    disp('volume ratio of leaf and macrogrid differing, please check!');
    res = 0;
  end;

end
  
