classdef feminfo < handle
%classdef feminfo < handle
%
% structure representing the fem-space information shared by all
% fem-functions. Implemneted as handle class, in order to be linked
% into df-classes.

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


% B. Haasdonk 13.1.2011

% adapted to new assembly

% I. Maier 24.03.2011

properties
   pdeg;
   dimrange;
   global_dof_index;
   ndofs;
   ndofs_per_element;

   % object of type gridbase test
   grid;
   nlagrange_nodes;
   lagrange_nodes;
   lagrange_edges;
   detDF;
   dirichlet_gids;
   robin_elind;
   robin_edgeind;
   neumann_elind;
   neumann_edgeind;
   l2_inner_product_matrix;
   h10_inner_product_matrix;
   regularized_h10_inner_product_matrix;
end

methods

  function df_info = feminfo (model,grid)
  %function feminfo(model,grid)
  % required fields of model:
  %       pdeg: polynomial degree
    df_info.pdeg = model.pdeg;
    df_info.dimrange = model.dimrange;
    df_info.ndofs = fem_ndofs(model,grid);
    df_info.ndofs_per_element = fem_ndofs_per_element(model);
    df_info.global_dof_index = fem_global_dof_index(model,grid);
    df_info.grid = grid;
    
    df_info.lagrange_nodes = lagrange_nodes_lcoord(model.pdeg);
    df_info.nlagrange_nodes = size(df_info.lagrange_nodes,1);
    df_info.lagrange_edges = lagrange_nodes_edges(model.pdeg);
    df_info.detDF = grid.A(:)*2; % determinant of reference mapping per element; 
    df_info.dirichlet_gids = [];

    % l2_inner_product and h10_inner product for now only for
    % scalar functions. Extend later.
    % make sure, that the folliwing is called with empty
    % dirichlet nodes, as fem_matrix_volume_component_assembly does
    % dirichlet treatment!!!
    % reaction term == 1 then reaction matrix is l2matrix
    if df_info.dimrange == 1
      model.reaction = @(grid,eindices,loc,params) ones(size(eindices,1),1);
      model.qdeg = 2*df_info.pdeg;
      df_info.l2_inner_product_matrix = ...
          fem_matrix_volume_part_assembly(@fem_matrix_reac_integral_kernel,model,df_info);
      %     H10_inner_product_matrix: identity diffusion
      model.diffusivity_tensor = @(grid,eindices,loc,params) ...
          [ones(size(eindices,1),1), zeros(size(eindices,1),1), ...
	   zeros(size(eindices,1),1), ones(size(eindices,1),1)]; 
      df_info.h10_inner_product_matrix = ...
          fem_matrix_volume_part_assembly(@fem_matrix_diff_integral_kernel,model,df_info);  
    end;
    
    % find dirichlet indices as preparation
    is_dirichlet_gid = zeros(df_info.ndofs,1);  
    % dirichlet-treatment: insert unit-vector row into A in dirichlet
    % rows. Little loop over local lagrange node index
    % this is repeated in r assembly as these will be kept separated
    for i=1:df_info.nlagrange_nodes
      % search all elements, in which the i-th Lagrange node is
      % dirichlet point
      eids = find(df_info.lagrange_edges(i,:)==1);
      if ~isempty(eids)
        for e = eids; % possibly one point belongs to 2 edges!!
          elids = find(grid.NBI(:,e)==-1);
          gids = df_info.global_dof_index(elids,i);
          is_dirichlet_gid(gids) = 1;
        end;
      end;
    end
    df_info.dirichlet_gids = find(is_dirichlet_gid);
    
    % find robin boundary elements and edge numbers
    ind = find(grid.NBI==-3); 
    [df_info.robin_elind,df_info.robin_edgeind] = ind2sub(size(grid.NBI),ind);
    
    % find neumann boundary elements and edge numbers
    ind = find(grid.NBI==-2); 
    [df_info.neumann_elind,df_info.neumann_edgeind] = ind2sub(size(grid.NBI),ind);
    % regularized_h10_inner_product_matrix: make unit rows in
    % dirichlet indices:
    
    if df_info.dimrange == 1
      reg_mat = df_info.h10_inner_product_matrix;
      if ~isempty(df_info.dirichlet_gids) 
        reg_mat(df_info.dirichlet_gids,:) = 0;
        A_dirichlet = sparse(df_info.dirichlet_gids,...
			     df_info.dirichlet_gids, ...
			     ones(size(df_info.dirichlet_gids)),...
			     df_info.ndofs,df_info.ndofs);
        reg_mat = reg_mat + A_dirichlet;
      end;
      df_info.regularized_h10_inner_product_matrix = reg_mat;
    end;
    
  end

end

end


