function v = fem_operators_output_boundary_integral(model,model_data)
%function v = fem_operators_output_boundary_integral(model,model_data)
%
% function computing the output vectors of an elliptic problem with
% finite element discretization and output functional consisting of
% a weighted boundary integral
%
%    `s(u) = \int_{\partial \Omega} f(x) u(x)`
%
% required fields of model:
%    output_function:  weight function `f`
%    decomp_mode:      decision about the mode to be performed 
%                      components, coefficients or complete

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

if model.decomp_mode == 2   % == 2 coefficients
  v = model.output_function([],model);
else
  % select all boundary elements for integration
  i = find(model_data.grid.NBI(:)<0);
  [elind,edgeind] =  ind2sub(size(model_data.grid.NBI),i);
  % turn off dirichlet deletion:
  with_dirichlet_deletion = 0;
  % assemble vectors
  v = fem_rhs_boundary_part_assembly(...
      @bnd_output_int_kernel,model,model_data.df_info,...
      elind,edgeind,with_dirichlet_deletion);
end;

%%%%%%%%%%%%%%%%%%%%% boundary integral kernel
function res = bnd_output_int_kernel(x,model,df_info,i,elind,edgeind)
% here x is a scalar value on unit interval
% f(x) = output_function * hatphi_i  
% multiplication with |edge| is realized in caller after quadrature.
% function can handle cell-array valued data
% this function is 99% identical to fem_rhs_neumann_integral_kernel
% perhaps merge sometime.
res = zeros(0,1);
for local_edge_id = 1:3
  lcoord = llocal2local(...
      df_info.grid,local_edge_id,x); % local coord for x on edge 1    
  hat_phi_i = fem_evaluate_basis_function(df_info,lcoord,i);
  inds = find(edgeind==local_edge_id);
  glob = local2global(df_info.grid,elind(inds),lcoord);
  w = model.output_function(glob,model);
  if ~iscell(w)
    res = [res;...
	   w * hat_phi_i; ...
	  ];
  else
    if ~iscell(res)
      res = cell(1,length(w));
      for q = 1:length(w)
	res{q} = zeros(0,1);
      end;
    end;
    for q = 1:length(w)
      res{q} = [res{q};...
		w{q} * hat_phi_i; ...
	       ];
    end;
  end;
end;

