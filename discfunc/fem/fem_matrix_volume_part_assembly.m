function A_comp = fem_matrix_volume_part_assembly(...
    A_int_kernel,model,df_info);
%function A_comp = fem_matrix_volume_part_assembly(...
%    A_int_kernel,model,df_info);
%
% auxiliary function assembling the volume integral components of
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

% B. Haasdonk 29.1.2011

A_comp = spalloc(df_info.ndofs,df_info.ndofs,1); 
for i = 1:df_info.nlagrange_nodes
  for j = 1:df_info.nlagrange_nodes
    gids_i = df_info.global_dof_index(1:df_info.grid.nelements,i);
    gids_j = df_info.global_dof_index(1:df_info.grid.nelements,j);
    f = @(x) A_int_kernel(x,model,df_info,i,j);      
    res = triaquadrature(model.qdeg,f); 
    if ~iscell(res)
      res = res .* df_info.detDF; % due to transformation formula
      A_tmp = sparse(gids_i,gids_j,res,df_info.ndofs,df_info.ndofs);
      A_comp = A_comp + A_tmp;
      %    if (A_tmp(6,11)~=0)
      %      A_tmp(6,11)
      %    end;
    else % iscell!!!
      if ~iscell(A_comp)	 
	% for first time: initialize A_comp
	A_comp = cell(1,length(res));
	for q=1:length(res)
	  A_comp{q} = spalloc(df_info.ndofs,df_info.ndofs,1); 
	end;
      end;
      for q = 1:length(res)
	res{q} = res{q} .* df_info.detDF; % due to transformation formula
	A_tmp = sparse(gids_i,gids_j,res{q},df_info.ndofs,df_info.ndofs);
	A_comp{q} = A_comp{q} + A_tmp;
      end;
    end;
  end;
end;  

if ~isempty(df_info.dirichlet_gids) 
  if ~iscell(A_comp)
    A_comp(df_info.dirichlet_gids,:) = 0;
  else
    for q = 1:length(A_comp)
      A_comp{q}(df_info.dirichlet_gids,:) = 0;
    end;
  end;
end;
