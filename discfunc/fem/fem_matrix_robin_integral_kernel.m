function res = fem_matrix_robin_integral_kernel(x,model,df_info,i,j,robin_elind,robin_edgeind)
%function res = fem_matrix_robin_integral_kernel(x,model,df_info,i,j,robin_elind,robin_edgeind)
%
% auxiliary function for integral kernel for A_robin
% integral kernel on all robin-edges simultaneously
% here x is a scalar value on unit interval
% f(x) = (alpha/beta + b^T n) hatphi_i hatphi_j 
% multiplication with |edge| is realized in caller after
% quadrature.
% This function can handle cell-array valued data

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
grid = df_info.grid;
for local_edge_id = 1:3
  lcoord = llocal2local(grid,local_edge_id,x); % local coord for x on edge 1    
  hat_phi_i = fem_evaluate_basis_function(df_info,lcoord,i);
  hat_phi_j = fem_evaluate_basis_function(df_info,lcoord,j);
  inds = find(robin_edgeind==local_edge_id);
  elinds = robin_elind(inds);
  b = model.velocity(grid,elinds,lcoord,model); % local model!

  if ~iscell(b)    
    res = [res; ( (model.robin_beta(grid,elinds,local_edge_id,x, ...
				    model).^(-1)) .* ...
		  model.robin_alpha(grid,elinds,local_edge_id,x, ...
				    model) ...
		  + b(:,1) .* grid.NX(elinds,local_edge_id) ...
		  + b(:,2) .* grid.NY(elinds,local_edge_id) ...
		  ) * hat_phi_i * hat_phi_j ]; % local model!
  else % iscell => components mode 
  
    % number of components = 1+ Q_velocity y 
    if ~iscell(res)
      res = cell(1,length(b)+1);
      for q = 1:length(b)+1
	res{q} = zeros(0,1);
      end;
    end;
    res{1} = [res{1}; ...
	      ( (model.robin_beta(grid,elinds,local_edge_id,x, ...
				  model).^(-1)) .* ...
		model.robin_alpha(grid,elinds,local_edge_id,x, ...
				  model) ...
		) * hat_phi_i * hat_phi_j ]; % local model!
    for q = 1:length(b)
      res{q+1} = [res{q+1};...
		  ( b{q}(:,1).*grid.NX(elinds,local_edge_id) ...
		    + b{q}(:,2).*grid.NY(elinds,local_edge_id) ...
		    ) * hat_phi_i * hat_phi_j; ...
		 ];
    end;
%    keyboard;  

  end;
end;
