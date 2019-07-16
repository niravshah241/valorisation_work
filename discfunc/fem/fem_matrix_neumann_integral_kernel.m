function res = fem_matrix_neumann_integral_kernel(x,model,df_info,i,j,neumann_elind,neumann_edgeind)
%function res = fem_matrix_neumann_integral_kernel(x,model,df_info,i,j,neumann_elind,neumann_edgeind)
%
% auxiliary function for integral kernel for A_neumann
% integral kernel on all neumann-edges simultaneously
% here x is a scalar value on unit interval
% f(x) = (b^T n) hatphi_i hatphi_j 
% multiplication with |edge| is realized in caller after quadrature.
% function can handle cell-array valued data

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

res = zeros(0,1);
for local_edge_id = 1:3
  % transform 1d coordinate to 2d local coordinate
  lcoord = llocal2local(df_info.grid,local_edge_id,x); % local coord for x on edge 1    
  hat_phi_i = fem_evaluate_basis_function(df_info,lcoord,i);
  hat_phi_j = fem_evaluate_basis_function(df_info,lcoord,j);
  inds = find(neumann_edgeind==local_edge_id);
  b = model.velocity(df_info.grid,neumann_elind(inds),lcoord,model); % local model!
  if ~iscell(b)
    res = [res;...
	   ( b(:,1).*df_info.grid.NX(neumann_elind(inds),local_edge_id) ...
	     + b(:,2).*df_info.grid.NY(neumann_elind(inds),local_edge_id) ...
	     ) * hat_phi_i * hat_phi_j; ...
	  ];
  else % iscell!!
    if ~iscell(res)
      res = cell(1,length(b));
      for q = 1:length(b)
	res{q} = zeros(0,1);
      end;
    end;
    for q = 1:length(b)      
      res{q} = [res{q};...
	     ( b{q}(:,1).*df_info.grid.NX(neumann_elind(inds),local_edge_id) ...
	     + b{q}(:,2).*df_info.grid.NY(neumann_elind(inds),local_edge_id) ...
	     ) * hat_phi_i * hat_phi_j; ...
	  ];
    end;
  
  end;
end;
