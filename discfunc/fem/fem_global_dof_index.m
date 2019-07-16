function global_dof_index = fem_global_dof_index(params,grid)
%function global_dof_index = fem_global_dof_index(params,grid)
% function computing the local-to-global dof map of a fem discrete function
%
% 'gid = global_dof_index(elid,lagrange_node)' yields the global
% index of the first dof of the basis function corresponding to the given
% Lagrange node and element 'elid'. 'gid:(gid+dimrange-1)' are the subsequent
% dof indices of the vectorial function in the lagrange node.  first all dofs
% in nodes are counted, then all dofs in element interior, then the dofs on
% edge-interiors.
%
% The Lagrange nodes 'l_1,...,l_m' with 'm=0.5*(pdeg+1)*(pdeg+2)' are sorted in
% the following order
% @verbatim
%       l_m = v_3
%       *
%       |\
%       | \
%       |  \
%       *   *
%       |    \
%       |     \
%       |______\
%       *   *   *
% v_1 = l_1      v_2 = l_(pdeg+1)
% @endverbatim
%
% where 'v_1,v_2,v_3' denote the sorting of the triangles corners.

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


% Bernard Haasdonk 11.1.2011


%if ~isfield(params,'debug')
%params.debug = 0;
%end;

    % generate global dof index map for scalar function

    nlagrange_nodes = (params.pdeg+1)*(params.pdeg+2)/2;
    global_dof_index_scalar = zeros(grid.nelements,nlagrange_nodes);   
    lind_finished = zeros(1,nlagrange_nodes);

    % 3 Lagrangenodes correspond to corners, 
    % global dof numbers: 1:nvertices
    lind_v1 = 1;    
    lind_v2 = 1+params.pdeg;    
    lind_v3 = (params.pdeg+1)*(params.pdeg+2)/2;   
    linds_corner = [lind_v1,lind_v2,lind_v3];
    global_dof_index_scalar(:,linds_corner) = ...
	grid.VI;
    lind_finished(linds_corner) = 1;
    
    % Lagrangenodes of interior: for each element
    % global dof numbers: (nvertices+1):(nvertices+nlangr_interior*nelements)
    lagrange_nodes = lagrange_nodes_lcoord(params.pdeg);
    linds_interior = find(...
	((1-lagrange_nodes(:,1)-lagrange_nodes(:,2))>eps) & ... 
	(lagrange_nodes(:,1)>eps) & ... 
	(lagrange_nodes(:,2)>eps) ... 
	);
    
    if ~isempty(linds_interior)
      gidtmp = (1:(length(linds_interior)*grid.nelements))+grid.nvertices;
      gidtmp = reshape(gidtmp,length(linds_interior),grid.nelements)';
      global_dof_index_scalar(:,linds_interior) = gidtmp;
      lind_finished(linds_interior) = 1;
    end;  
    
    % Lagrangenodes on edges: are remaining ones
    % global dof numbers: 
    %       nvertices + nlangr_interior*nelements + 1
    % until nvertices + nlangr_interior*nelements + nedges_interior * (pdeg-1)
    linds_edges = cell(3,1);
    linds_edges{1} = find( ...
    	(lagrange_nodes(:,2)<eps) & ... 
    	((1-lagrange_nodes(:,1)-lagrange_nodes(:,2))>eps) & ... 
    	(lagrange_nodes(:,1)>eps) ... 
      );

    linds_edges{2} = find( ...
    	((1-lagrange_nodes(:,1)-lagrange_nodes(:,2))<eps) & ... 
    	(lagrange_nodes(:,1)>eps) & ... 
    	(lagrange_nodes(:,2)>eps) ... 
	);

    linds_edges{3} = find( ...
    	(lagrange_nodes(:,1)<eps) & ... 
    	((1-lagrange_nodes(:,1)-lagrange_nodes(:,2))>eps) & ... 
    	(lagrange_nodes(:,2)>eps) ... 
	);
    % switch order of 3rd edge dofs, such that all are counterclockwise
    linds_edges{3} = linds_edges{3}(end:-1:1);
    
    % insert new nodes for edges with element index larger than
    % neighbour => boundary and inner edges treated once.    
    elids = cell(3,1);
    offset = (length(linds_interior)*grid.nelements)+grid.nvertices;
    for i=1:3
      elids{i} = find(grid.NBI(:,i)<(1:grid.nelements)');
      gidtmp = (1:(length(linds_edges{i})*length(elids{i})))+offset;
      gidtmp = reshape(gidtmp,length(linds_edges{i}),length(elids{i}))';
      global_dof_index_scalar(elids{i},linds_edges{i}) = gidtmp;
      offset = offset + (length(linds_edges{i})*length(elids{i}));    
    end;
    
    % set correct gid for elements with index smaller than neighbour
    not_elids = cell(3,1);
    lind_inv_edges = cell(3,1); % invert order
    for i = 1:3
      % not_elids{i} satisfies
      % grid.NBI(not_elids{i},i)>not_elids{i});
      % hence, not_elids{i} are the number of the
      % elements, with index smaller than neighbour
      % and the nieghbours are accessable over local edge number i
      not_elids{i} = find(grid.NBI(:,i)>(1:grid.nelements)');
      lind_inv_edges{i} = linds_edges{i}(end:-1:1);
    end;

    % loop over all different local edge number combinations
    for j = 1:3 % edge number on lower element index side   
      not_elids_neighbour = grid.NBI(not_elids{j},j);
%      if params.debug
%	if ~isequal(grid.NBI(not_elids_neighbour,j),not_elids{j});
%	  disp('problem in dof enumeration!');
%	  keyboard;
%	end;
%      end;
      for i = 1:3 % edge number on larger element index side
	subind_i = find(grid.NBI(not_elids_neighbour,i)==not_elids{j});
        gidtmp = global_dof_index_scalar(...
	    not_elids_neighbour(subind_i), ...
	    linds_edges{i});
            global_dof_index_scalar(...
		not_elids{j}(subind_i),...
		lind_inv_edges{j}) = ...
		gidtmp;
      end;
    end;
    
    % mark finished
    for i = 1:3
      lind_finished(linds_edges{i}) = 1;
    end;
        
    % checks:
    if ~isempty(find(lind_finished==0))
      error('not all lagrange indices treated!!!');
    end;
    
    if ~isempty(find(global_dof_index_scalar == 0))
      error('not all global dof numbers set!!');
    end;
    
    % map for vectorial function simple scaling keeping 1 fix:
    global_dof_index = (global_dof_index_scalar -1)* params.dimrange+1;

    % check: local dofs mapped to identical global dofs must have
    % identical global coordinates of lagrange point.    
    %disp('implement test!');
    %keyboard;
    
