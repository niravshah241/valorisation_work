function [A_diff , A_adv, A_reac, A_dirichlet, A_neumann, A_robin] = ...
    fem_matrix_parts_assembly(model,df_info)
%function [A_diff , A_adv, A_reac, A_dirichlet, A_neumann, A_robin] = ...
%    fem_matrix_parts_assembly(model,df_info)
%
% function assembling the system matrix parts for solving an
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

% Bernard Haasdonk 13.1.2011

% modification: data functions of the 'model' are now assumed to be local
% functions, see elliptic_discrete_model

% I. Maier 24.03.2011

% assemble system
% A = (a_ij)_i,j=1..ndofs
%
% a_i. = unit_row_vector(i)                         if x_i is in Gamma_D
%
% a_ij = a_ij_diff + a_ij_adv + a_ij reac 
%      + a_ij_robin + a_ij_neumann                  otherwise
%
% 
% a_ij_diff = \int_Omega (\grad phi_j)^T * a(x)^T * \grad \phi_i 
%
% a_ij_adv = - \int_Omega \phi_j * (b^T \div \phi_i)
%         caution: minus sign!!!!
%
% a_ij_reac = \int_Omega c phi_j phi_i
%
% a_ij_neumann = \int_Gamma_N (b^T n) * phi_j * phi_i 
%
% a_ij_robin =  \int_Gamma_R (alpha/beta + b^T*n) * phi_j * phi_i

% A_diff
% a_ij = int_Omega (grad phi_j)^T * A^T * \grad phi_i
% by simultaneous integration over all elements and
% transformation  formula for reference triangle, i.e. 
% multiplication with detDF and Jacobianinversetransposed (JIT)
%
% int_T (grad phi_j)^T * A^T * grad phi_i
% = int_hatT  (\grad hatphi_l)^T * JIT^T * A^T * (JIT * grad hatphi_k) detDF 
% with suitable local indices k and l
%
% function supports affine decomposition, i.e. depending on model.decomp_mode
% This program is open source.  For license terms, see the COPYING file.

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


if model.decomp_mode == 2 % == coefficients
 
  A_diff = 0;
  if model.has_diffusivity
    A_diff = model.diffusivity_tensor([],[],[],model);
  end;
  
  vel =[];
  A_adv = 0;
  if model.has_advection
    % minus sign has been respected in components, no need here!
    vel = model.velocity([],[],[],model);
    A_adv = vel;
  end;

  A_robin = 0;
  if model.has_robin_values
    A_robin = [1];
    if model.has_advection
      A_robin = [1,vel(:)];
    end;
  end;
  
  A_neumann = 0;
  if (model.has_neumann_values) & (~isempty(vel));
    A_neumann = vel;
  end;
  
  A_reac = 0;
  if model.has_reaction
    A_reac = model.reaction([],[],[],model);
  end;

  A_dirichlet = 1; % always one component, coefficient 1
    
else % decomp_mode == 0 (complete) or 1 (components)
  
  ndofs = df_info.ndofs;
  nelements = df_info.grid.nelements;

  zero_mat = spalloc(ndofs,ndofs,1);
  if model.decomp_mode == 1
    zero_mat = {zero_mat};
  end;
  
  A_diff =zero_mat;
  if model.has_diffusivity
    A_diff = fem_matrix_volume_part_assembly(...
	@fem_matrix_diff_integral_kernel,model,df_info);
  end;
  
  % A_adv: !!!!!!!!!!!!! note minus sign here !!!!!!!!!!!!
  A_adv = zero_mat;
  if model.has_advection
    A_adv = fem_matrix_volume_part_assembly(...
	@fem_matrix_adv_integral_kernel,model,df_info);
    if ~iscell(A_adv)
      A_adv = -A_adv;
    else
        for q = 1:length(A_adv)
          A_adv{q} = - A_adv{q};
        end;    
    end;
  end;
  
  % A reac
  A_reac = zero_mat;
  if model.has_reaction
    A_reac = fem_matrix_volume_part_assembly(...
	@fem_matrix_reac_integral_kernel,model,df_info);
  end;
  
  % A_neumann
  A_neumann = zero_mat;
  if model.has_advection
    A_neumann = fem_matrix_boundary_part_assembly(...
	@fem_matrix_neumann_integral_kernel,model,df_info,df_info.neumann_elind,df_info.neumann_edgeind);
  end;
  
  % A_robin
  A_robin = zero_mat;
  if model.has_robin_values
    A_robin = fem_matrix_boundary_part_assembly(...
	@fem_matrix_robin_integral_kernel,model,df_info,df_info.robin_elind,df_info.robin_edgeind);
  end;
  
  % A_dirichlet
  if ~isempty(df_info.dirichlet_gids)
    A_dirichlet = sparse(df_info.dirichlet_gids,...
			 df_info.dirichlet_gids, ...
			 ones(size(df_info.dirichlet_gids)),...
			 ndofs,ndofs);
  else
    A_dirichlet = spalloc(ndofs,ndofs,1);
  end;
 
  if model.decomp_mode == 1
    A_dirichlet = {A_dirichlet}; % make cell array.
  end;
  
end;
  
	  
