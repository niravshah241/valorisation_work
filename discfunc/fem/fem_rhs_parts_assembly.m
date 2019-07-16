function [r_source, r_dirichlet, r_neumann, r_robin] = ...
    fem_rhs_parts_assembly(model,df_info)
%function [r_source, r_dirichlet, r_neumann, r_robin] = ...
%    fem_rhs_parts_assembly(model,df_info)
%
% function assembling the right hand side vector parts for solving an
% elliptic pde with finite elements
%
% - div ( a(x) grad u(x)) + div (b(x) u(x)) + c(x) u(x) = f(x)    on Omega
%                                                 u(x)) = g_D(x)  on Gamma_D
%                                   a(x) (grad u(x)) n) = g_N(x)  on Gamma_N
%          alpha(x) u(x) +  beta(x) a(x) (grad u(x)) n) = g_R(x)  on Gamma_R
%
%                                   s = l(u) linear output functional
% 
%  Here, we denote the functions as
%                   u: scalar 'solution' (if known, for validation purpose)
%                   f: scalar 'source'
%                   a: tensor valued 'diffusivity_tensor'
%                   b: vector valued 'velocity'
%                   c: scalar 'reaction'
%                 g_D: scalar 'dirichlet_values'
%                 g_N: scalar 'neumann_values'
%                 g_R: scalar 'robin_values'
%               alpha: scalar 'robin_alpha'
%                beta: scalar 'robin_beta'
%
% which are contained in a 'model'.See poisson_model for an
% extensive description. 
% df_info is discrete-function-info containing the grid
% global_dof_index, dirichlet_indices, lagrange-nodes and further
% fem-information used by the system and the rhs assembly by
% fem_rhs_parts_assembly
%
% function allows affine decomposition, i.e. result depending on 
% model.decomp_mode

% Bernard Haasdonk 13.1.2011

% modification: data functions of the 'model' are now assumed to be local
% functions, see elliptic_discrete_model

% I. Maier 24.03.2011


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

if model.decomp_mode == 2

  r_source = 0;
  if model.has_source
    r_source = model.source([],[],[],model);
  end;

  r_dirichlet = 0;
  if model.has_dirichlet_values
    r_dirichlet = model.dirichlet_values([],[],[],[],model);
  end;

  r_neumann = 0;
  if model.has_neumann_values
    r_neumann = model.neumann_values([],[],[],[],model);
  end;

  r_robin = 0;
  if model.has_robin_values
    r_robin = model.robin_values([],[],[],[],model);
  end;

else % decomp_mode == 1 or 0
    
  ndofs = df_info.ndofs;
  nelements = df_info.grid.nelements;
  nlagrange_nodes = df_info.nlagrange_nodes;
  grid = df_info.grid;
  qdeg = model.qdeg;
  llcoord = lagrange_nodes_edges_llcoord(model.pdeg); % 1d coordinates of lagrange nodes

  if model.decomp_mode == 0
    zerovec = zeros(ndofs,1);
  else
    zerovec = {zeros(ndofs,1)};
  end;
  
  % right hand side vector assembly
  %
  % r = (r_i)_i=1...ndofs
  % 
  % r_i = g_D(x_i)    if x_i in \Gamma_D
  %
  % r_i = r_i_source + r_i_neumann + r_i_robin
  %
  %       r_i_volume = \int_Omega q * phi_i
  %
  %       r_i_neumann = \int_Gamma_N g_N * phi_i
  %
  %       r_i_robin = \int_Gamma_R g_R / beta * phi_i
  
  % volume integral  
  r_source = zerovec;
  if model.has_source
    r_source = fem_rhs_volume_part_assembly(...
	model.source,model,df_info);
  end;
  
  % neumann terms
  r_neumann = zerovec;
  if model.has_neumann_values  
    r_neumann = fem_rhs_boundary_part_assembly(...
	@fem_rhs_neumann_integral_kernel,model,df_info,df_info.neumann_elind,df_info.neumann_edgeind);
  end;
  
  % robin terms
  r_robin = zerovec;
  if model.has_robin_values  
    r_robin = fem_rhs_boundary_part_assembly(...
	@fem_rhs_robin_integral_kernel,model,df_info,df_info.robin_elind,df_info.robin_edgeind);
  end;
  
  r_dirichlet = zerovec; 
  if model.has_dirichlet_values
    % dirichlet-treatment: insert dirichlet values into rhs
    % little loop over local lagrange node index
    % this will be repeated in A assembly as these will be kept separated
    for i=1:nlagrange_nodes
      % search all elements, in which the i-th Lagrange node is
      % dirichlet point
      eids = find(df_info.lagrange_edges(i,:)==1);
      if ~isempty(eids)
	for e = eids; % possibly one point belongs to 2 edges!!
	  elids = find(grid.NBI(:,e)==-1);
	  %if llcoord(i,e) == -1
	  %  disp('wrong llcoord values in lagrange_nodes_edges_llcoord!');
	  %end;
	  dirvals = model.dirichlet_values(grid,elids,e,llcoord(i,e),model); % local model!
	  gids = df_info.global_dof_index(elids,i);
	  if model.decomp_mode == 0
	    r_dirichlet(gids) = dirvals;
	  else % decomp_mode = 1
	    if length(dirvals)~=length(r_dirichlet)
	      r_dirichlet = cell(1,length(r_dirichlet));
	      for q = 1:length(dirvals)
		r_dirichlet{q} = zeros(ndofs,1);
	      end;
	    end;
	    for q = 1:length(dirvals)
	      r_dirichlet{q}(gids) = dirvals{q};
	    end;	    
	  end;
	end;
      end;
    end;  
  end;
 
end;

  



