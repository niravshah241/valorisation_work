function r_comp = fem_rhs_boundary_part_assembly(...
    r_int_kernel,model,df_info,elind,edgeind,with_dirichlet_deletion);
%function r_comp = fem_rhs_boundary_part_assembly(...
%    r_int_kernel,model,df_info,elind,edgeind,with_dirichlet_deletion);
%
% auxiliary function assembling the boundary integral components of
% system matrix A, i.e. Neumann and Robin components
%
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

r_comp = zeros(df_info.ndofs,1);
for i = 1:df_info.nlagrange_nodes
  EL = [];
  for local_edge_ind = 1:3;
    ind = find(edgeind==local_edge_ind);
    EL = [EL; df_info.grid.EL(elind(ind),local_edge_ind)];
  end;
  gids = df_info.global_dof_index(elind,i);
  f = @(x) r_int_kernel(x,model,df_info,i,elind,edgeind);
  res = intervalquadrature(model.qdeg,f);

  if ~iscell(res)
    res = res .* EL; % due to transformation formula
    
    % here sparse matrix must be used to correclty handle multiple
    % identical gids accumulatively!!!
    r_comp_inc = sparse(gids,ones(length(gids),1),res,df_info.ndofs,1);
    r_comp = r_comp + r_comp_inc;
    
  else % iscell!!
    if ~iscell(r_comp)	 
      % for first time: initialize r_comp
      r_comp = cell(1,length(res));
      for q=1:length(res)
	r_comp{q} = zeros(df_info.ndofs,1); 
      end;
    end;
    for q = 1:length(res)
      res{q} = res{q} .* EL; % due to transformation formula
      r_comp_inc = sparse(gids,ones(length(gids),1),res{q},df_info.ndofs,1);
      r_comp{q} = r_comp{q} + r_comp_inc;
    end; 
  end;
end;  

if nargin<6
  with_dirichlet_deletion = 1;
end;

if with_dirichlet_deletion
  if ~isempty(df_info.dirichlet_gids)
    if ~iscell(r_comp)
      r_comp(df_info.dirichlet_gids) = 0;
    else
      for q = 1:length(r_comp)
	r_comp{q}(df_info.dirichlet_gids) = 0; 
      end;
    end;
  end;
end;
