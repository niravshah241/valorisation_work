function res = fem_rhs_robin_integral_kernel(x,model,df_info,i,robin_elind,robin_edgeind);
%function res = fem_rhs_robin_integral_kernel(x,model,df_info,i,robin_elind,robin_edgeind);
%
% auxiliary function for integral kernel for r_robin
% integral kernel on all robin-edges simultaneously
% here x is a scalar value on unit interval
% f(x) = g_R/beta * hatphi_i  
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
  lcoord = llocal2local(df_info.grid,local_edge_id,x); % local coord for x on edge 1    
  hat_phi_i = fem_evaluate_basis_function(df_info,lcoord,i);
  inds = find(robin_edgeind==local_edge_id);
  beta = model.robin_beta(df_info.grid,robin_elind(inds),local_edge_id,x,model); % local model!
  g_R = model.robin_values(df_info.grid,robin_elind(inds),local_edge_id,x,model); % local model!
  if ~iscell(g_R)
    res = [res;...
	   (beta.^(-1)).*g_R * hat_phi_i; ...
	];
  else % iscell!!
    if ~iscell(res)
      res = cell(1,length(g_R));
      for q = 1:length(g_R);
	res{q} = zeros(0,1);
      end;
    end;
    for q = 1:length(g_R)
      res{q} = [res{q};...
		(beta.^(-1)).*g_R{q} * hat_phi_i; ...
	       ];
    end;
  end;
end;

