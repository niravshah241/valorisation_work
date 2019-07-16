function A_comp = fem_matrix_boundary_part_assembly(...
    A_int_kernel,model,df_info,elind,edgeind);
%function A_comp = fem_matrix_boundary_part_assembly(...
%    A_int_kernel,model,df_info,elind,edgeind);
% auxiliary function assembling the boundary integral components of
% system matrix A
% note: cell-array valued kernels can be integrated.

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

% B. Haasdonk 22.2.2011

A_comp = spalloc(df_info.ndofs,df_info.ndofs,1); 
for i = 1:df_info.nlagrange_nodes
  for j = 1:df_info.nlagrange_nodes
    f = @(x) A_int_kernel(x,model,df_info,i,j,elind,edgeind);
    res = intervalquadrature(model.qdeg,f);
    EL = [];
    for local_edge_ind = 1:3;
      ind = find(edgeind==local_edge_ind);
      EL = [EL; df_info.grid.EL(elind(ind),local_edge_ind)];
    end;
    gids_i = df_info.global_dof_index(elind,i);
    gids_j = df_info.global_dof_index(elind,j);
    
    if ~iscell(res)
      res = res .* EL; % due to transformation formula
      A_tmp = sparse(gids_i,gids_j,res,df_info.ndofs,df_info.ndofs);
      A_comp = A_comp + A_tmp;
    else % iscell!!!      
      if ~iscell(A_comp)	 
	% for first time: initialize A_comp
	A_comp = cell(1,length(res));
	for q=1:length(res)
	  A_comp{q} = spalloc(df_info.ndofs,df_info.ndofs,1); 
	end;
      end;
      for q = 1:length(res)
	res{q} = res{q} .* EL; % due to transformation formula
	A_tmp = sparse(gids_i,gids_j,res{q},df_info.ndofs,df_info.ndofs);
	A_comp{q} = A_comp{q} + A_tmp;
      end;      
    end;
  end;
end;

if ~isempty(df_info.dirichlet_gids)  
  if ~iscell(res)
    A_comp(df_info.dirichlet_gids,:) = 0;
  else % iscell!
    for q = 1:length(A_comp)
      A_comp{q}(df_info.dirichlet_gids,:) = 0;
    end;
  end;
end;
