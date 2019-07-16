function ngrid = refine(grid, lids)
%function ngrid = refine(grid, lids)
% refinement routine for cubegrid.
%
% Refinement rule is a simple division into `2^dim` smaller cubes After
% refinement, a condensation of the vertex list is performed by
% remove_duplicate_vetrices(), i.e. duplicate vertices generated during
% refinement are removed.
%
% Parameters:
% lids: List of leaf-element ids of the to-be-refined codim '0' entities.  Only
%       leaf-elements can reasonably be refined, so do not pass global elements
%       here, as internally a translation into global element ids is performed.
%
% Return values:
%  ngrid: the refined grid

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
  
%  Make sure, that only leaf-elements are refined
  
  %  m = min(grid.isleaf(ids));
  %  if (m==0)
  %    error('Refinement of non-leaf entity is requested!');
  %  end;
    
  dim = grid.dimension;
  
  % get global ids
  gids = lid2gid(grid,lids);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % set up refinement rule for elements:
  % for each dimension, a weight matrix is set up, which is later
  % multiplied to vertex coordinates
  % 1D: 3 child points:     ref_weights = [1.0,    0.5,   0.0 ; ... 
  %                                        0.0,    0.5    1.0];
  %     2 new elements:     ref_indices = [1, 2;     
  %                                        2, 3];
  % 2D  9 child points :    ref_weights = 
  %         [1.0,    0.5   0.0     0.5  0.25 0.5  0.0   0.0  0.0; ... 
  %          0.0,    0.5   1..0    0.0  0.25 0.5  0.0   0.0  0.0; ...
  %          0.0,    0.0,  0.0     0.5  0.25 0.5  1.0   0.5  0.0; ...
  %          0.0,    0.0,  0.0     0.0  0.25 0.5  0.0   0.5  1.0]
  %     4 new elements (doppelt so viele wie bei 1D weniger) 
  %              ref_indices = [1 2 4 5 ; ...
  %                             2 3 5 6 ; ...
  %                             4 5 7 8 ; ...
  %                             5 6 8 9 ];

  % Rekursiver Aufbau der weights:   
  % OD: ref_weights[0] = [1.0]
  % 1D: ref_weights[1] = [ref_weights[0] ,   0.5*rw[0] ;    0.0...
  %                              0.0         0.5*rw[0] ,   rw[0] ]
  % analog weiter :-)
  
  rw = ones(1,1);
  for i=1:dim;
    rw = [rw,                0.5*rw,  zeros(size(rw)); ...
	  zeros(size(rw)),   0.5*rw,       rw        ];
  end;
  ref_weights = rw;
  
  % rekursiver Aufbau der indizes:
  %     shift = 3^(dim-1)
  % ri = [ ri, ri+ shift ; ...
  %	   ri+shift, ri + 2*shift]; 
  %
  ri = ones(1,1);
  for i = 1:dim;
    shift = 3^(i-1);
    ri = [ri, ri + shift; ...
	  ri+shift, ri+2*shift];
  end;
  ref_indices = ri;
  
  % generate new points after refinement:  
  % all 3^dim points in all to-be-refined elements
  
  new_vertex = zeros(0,dim);
  % generate new elements after refinement
  new_vertexindex = zeros(0,2^dim);
  new_level = [];
  new_firstchild = [];
  gids_firstchild = grid.nelements(1)+2^dim*(0:length(gids)-1)+1;
  
  for i = gids';
    co = grid.vertex(grid.vertexindex(i,:),:); % glob vertex coords as rows
    vid_offset = size(new_vertex,1) + size(grid.vertex,1);
    new_vertex = [new_vertex; ref_weights' * co]; 
    new_firstchild = [new_firstchild; zeros(2^dim,1)];
    new_vertexindex = [new_vertexindex; ref_indices + vid_offset];    
    new_level = [new_level; ones(2^dim,1)*grid.level(i)+ 1];
  end;
  
  % now consistent list is obtained, but vertices might be duplicated.
  grid.vertex = [grid.vertex; new_vertex];
  grid.vertexindex = [grid.vertexindex; new_vertexindex];
  grid.firstchild = [grid.firstchild; new_firstchild];
  grid.firstchild(gids) = gids_firstchild;
  grid.nelements = size(grid.vertexindex,1);
  grid.nvertices = size(grid.vertex,1); 
  grid.level = [grid.level(:); new_level];
  grid.isleaf = [grid.isleaf(:); ones(length(new_level),1)];
  % set refined elements to non-leaf
  grid.isleaf(gids) = 0; 
  
  % increase ref-counter 
  grid.refine_steps = grid.refine_steps + 1; 
  
  grid.creation_step = [grid.creation_step(:); ...
		    grid.refine_steps*ones(length(new_level),1)];
    
  ngrid = remove_duplicate_vertices(grid,1e-20);  
%  ngrid = grid;

  if ~check_consistency(ngrid)
    disp('problem in duplicate elimination, keeping duplicate vertices!')
    ngrid = grid;
  end;

end

