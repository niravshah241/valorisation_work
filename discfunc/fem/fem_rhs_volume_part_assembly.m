function r_comp = fem_rhs_volume_part_assembly(...
    r_int_kernel,model,df_info);
%function r_comp = fem_rhs_volume_part_assembly(...
%    r_int_kernel,model,df_info);
%
% auxiliary function assembling the volume integral components of
% right hand side, i.e. source components
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
% I. Maier 24.03.2011

r_comp = zeros(df_info.ndofs,1);

for i = 1:df_info.nlagrange_nodes
  basis_func_i = @(x) fem_evaluate_basis_function(df_info,x,i);
  gids = df_info.global_dof_index(1:df_info.grid.nelements,i);
  % integral kernel: q*phi on all elements simultaneously
  % and transformation on reference triangle, i.e. mult with detDF
  % here x is in reference triangle! auxiliary function defined below
  f = @(x) int_kernel_mult_phi_i(x,r_int_kernel,basis_func_i,model,df_info);
  
  res = triaquadrature(model.qdeg,f); 
  if ~iscell(res)
    res = res.* df_info.detDF; % due to transformation formula
    
    % caution: if multiple gids identical, only the last entry is
    % added in case of real vectors! 
    % => r_source(gids) = r_source(gids) + res_unique; !!!
    % Instead here sparse matrix is
    % used, there assignment is incremental
    % add all increments of identical gids!!!
    r_comp_incr = sparse(gids,ones(length(gids),1),res,df_info.ndofs,1);
    r_comp = r_comp + r_comp_incr;
    %    r_source_nterms_incr = sparse(gids,ones(length(gids),1),...
    %				  ones(length(gids),1),ndofs,1);
    %    r_source_nterms = r_source_nterms + r_source_nterms_incr;
  
  else % iscell!!!
    if ~iscell(r_comp)	 
      % for first time: initialize r_comp
      r_comp = cell(1,length(res));
      for q=1:length(res)
	r_comp{q} = zeros(df_info.ndofs,1); 
      end;
    end;
    for q = 1:length(res)
      res{q} = res{q}.* df_info.detDF; % due to transformation formula
      r_comp_incr = sparse(gids,ones(length(gids),1),res{q},df_info.ndofs,1);
      r_comp{q} = r_comp{q} + r_comp_incr;
    end;
  end;
end;

if ~isempty(df_info.dirichlet_gids)
  if ~iscell(r_comp)
    r_comp(df_info.dirichlet_gids) = 0;
  else
    for q = 1:length(r_comp)
      r_comp{q}(df_info.dirichlet_gids) = 0; 
    end;
  end;
end;


%%%%%%%%%%%%%%%%%%%%%%%%%% auxiliary function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function res = int_kernel_mult_phi_i(x,r_int_kernel,basis_func_i,model,df_info)

r_int = r_int_kernel(df_info.grid,1:df_info.grid.nelements,x,model);
if ~iscell(r_int)
  res = r_int * ...
	basis_func_i(x);
else
  res = cell(1,length(r_int));
  for q = 1:length(r_int)
    res{q} = r_int{q} * ...
	     basis_func_i(x);
    
  end;
end;
