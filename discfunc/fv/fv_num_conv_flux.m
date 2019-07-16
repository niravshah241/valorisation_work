function num_flux = fv_num_conv_flux(model,model_data,U)
%function num_flux = fv_num_conv_flux(model,model_data,U)
%
% function computing a numerical convective flux
% matrix. Dirichlet-boundary treatment is performed, but Neuman-values
% are set to nan, must be performed by calling function. Currently implemented
% lax-friedrichs and engquist-osher.
%
% fields of num_flux:
%   Lg: lipschitz constant satisfying 
%       |g_ij(w,v) - g_ij(w',v')| <= Lg |S_ij| ( |w'-w| + |v'-v| ) 
%   G:  the matrix of numerical flux values across the edges.
%       boundary treatment is performed as follows:  
%       in dirichlet-boundary points, the neighbour-values are simply set
%       to this value and the flux computed. 
%       Neuman-treatment is not performed here
%
% required fields of model:
%   name_convective_num_flux      : 'none', 'lax-friedrichs', 'engquist-osher'
%
% in case of lax-friedrichs, either as diffusivity the value of the
% flux_matrix is used or explicit setting of the value by
% model.lxf_lambda is possible.
%  
% for engquist-osher, we assume for the convective flux function:  
%         1. f'(u)*n does not change sign on interval(u,0) 
%         2. f(0) = 0
% additionally, a conv_flux_derivative data-function must be available.
%  
% plus additional fields required by dirichlet_values
%
% if grid is empty, it is generated

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

  
% Bernard Haasdonk 10.4.2006

  grid = model_data.grid;

%  if ~isfield(model,'flux_quad_degree') 
%    model.flux_quad_degree = 1;
%  end;
  
  num_flux = [];
  num_flux.G = nan* ones(size(grid.NBI));
  
  % determine index sets of boundary types
  real_NB_ind = find(grid.NBI>0);
  NBIind = grid.NBI(real_NB_ind);
  INBind = grid.INB(real_NB_ind);
  % NB_real_ind comprises comprises some kind of inverse cell-neighbour
  % relation that leads back to the original cell from where we started
  % Falls 
  %  [i,j]=ind2sub(size(grid.NBI),NB_real_NB_ind(sub2ind(size(grid.NBI),k,l)))
  % dann ist
  %  grid.NBI(i,j) == k
  NB_real_NB_ind = sub2ind(size(grid.NBI),NBIind,INBind);
  
  dir_NB_ind = find(grid.NBI == -1);
  
  % matrix with the neighbouring function values  
  % UNB = nan * ones(size(flux_mat.Fx));
  UNB = nan * ones(size(grid.NBI));
  UNB(real_NB_ind) = U(grid.NBI(real_NB_ind));  
  UU = repmat(U,1,size(grid.NBI,2));
    
  %%%%%%%%%%%%%%%%%%%  convective numerical flux: Lax-Friedrichs %%%%%%%%%%%%%%
  if isequal(model.name_convective_num_flux,'nonlinear-diffusion')   
    %  g_jl(u,v) = |S_ij| K(grid.ESX,grid.ESY, u(ES))(psi(v)-psi(u))
    %  z.B. u(ES) = 0.5(u+v)
    
    
  %%%%%%%%%%%%%%%%%%%  convective numerical flux: Lax-Friedrichs %%%%%%%%%%%%%%
  elseif isequal(model.name_convective_num_flux,'lax-friedrichs')   
    % Lax-Friedrichs-Flux in 2D:
    %   g_jl(u,v) = 1/2 * nu_jl (f(u,Sx,Sy) + f(v,Sx,Sy)) - 
    %             |S_jl|/(2 lambda_jl) (v-u)  
    %   or equivalently
    %   g_jl(u,v) = 1/2 |S_ij| * (n_jl * (f(u)+f(v)) - ...
    %                           1/(2 lambda_jl) (v-u))
    %   mit lambda_jl := dt / |m_j - m_l|
    
    % FxNB and FyNB and UNB are computed, where possible:
    %    both are valid in real_NB_ind
    
    % matrix evaluating the flux by quadratures over all edges
    flux_mat = fv_conv_flux_matrix(model,model_data,U);
    
    % in case of divergence cleaning, the velocity field must be
    % available!! This is obtianed by setting u to 1  
    if (model.divclean_mode>0)
      Udummy = ones(size(U));
      flux_mat_Uone = fv_conv_flux_matrix(model,model_data,Udummy);    
    end;
    
    % matrix with flux values, the neighbour values U(j) inserted
    % values must be reflected across all edges
    % i.e. GNB(i,j) = flux.G(NBI(i,j),INB(i,j))
    % take care of boundary-indices!!
    FxNB = nan * ones(size(flux_mat.Fx));
    %    [i.j] = ind2sub(size(flux_mat.G),real_NB_ind);
    FxNB(real_NB_ind) = flux_mat.Fx(NB_real_NB_ind);
    
    FyNB = nan * ones(size(flux_mat.Fy));
    %    [i.j] = ind2sub(size(flux_mat.G),real_NB_ind);
    FyNB(real_NB_ind) = flux_mat.Fy(NB_real_NB_ind);
    
    % L1 is the supremum in space of the derivative of the first component  
    % L2 is the supremum in space of the derivative of the second component  
    L1 = max(max(abs(flux_mat.Fx)));
    L2 = max(max(abs(flux_mat.Fy)));
  
    % Lipschitz-constant of LxF Flux
    num_flux.Lg = 1/(2*flux_mat.lambda) + 1/2 * norm([L1, L2]);
    % an alternative formula without motivation:
    %    num_flux.Lg =  max(max(grid.S)) * model.c * 0.25; 
  
    % in the linear case the following worked:
    % numflux_mat.G = 1/2 * grid.S .* ...
    %	        (grid.Nx.* (flux_mat.Vx.*(OU + OUext(grid.NBI))) + ...
    %	      grid.Ny.* (flux_mat.Vy.*(OU + OUext(grid.NBI)))) ...
    %	+ 1/(2 * flux_mat.lambda) * grid.S .* (OU - OUext(grid.NBI));     

    lambda = flux_mat.lambda;
    if ~isempty(model.lxf_lambda)
      lambda = model.lxf_lambda;
    end;
  
    % in the general nonlinear case:
    num_flux.G = 1/2 * grid.EL .* ...
	(grid.NX.* (flux_mat.Fx + FxNB) + ...
	 grid.NY.* (flux_mat.Fy + FyNB)) ...
    	+ 1/(2 * lambda) * grid.EL .* (UU - UNB);     

    % Dirichlet boundary treatment
    if ~isempty(dir_NB_ind>0)
    % determine dirichlet-boundary values as required by convective 
    % and diffusive flux. Neumann determined at end
      Xdir = grid.ECX(dir_NB_ind);
      Ydir = grid.ECY(dir_NB_ind);
      Udir = model.dirichlet_values_ptr([Xdir,Ydir],model);
      
      % if conv_flux_matrix has been divergence cleaned, no further
      % flux evaluations may be performed!!
      if ~isequal(model.divclean_mode,'none')
	if model.verbose>9
	  disp('skipping boundary flux eval. due to divergence cleaning!!');
	end;
	% the following is unstable, if U reaches 0!!!!
	%	Fdir.Fx = flux_mat.Fx(dir_NB_ind)./UU(dir_NB_ind).*Udir;
	%	Fdir.Fy = flux_mat.Fy(dir_NB_ind)./UU(dir_NB_ind).*Udir;
	% instead again evaluation of the whole matrix.
	Fdir.Fx = flux_mat_Uone.Fx(dir_NB_ind).*Udir;
	Fdir.Fy = flux_mat_Uone.Fy(dir_NB_ind).*Udir;
	
      else % ordinary flux evaluation by quadratures
	   %	Fdir = conv_flux(model,Udir, Xdir, Ydir);
	   [elids, edgeids] = ind2sub(size(grid.VI),dir_NB_ind);
	   PP = edge_quad_points(grid,elids, edgeids, ...
				 model.flux_quad_degree);
	   
	   % the following is only relevant in case of use of a 
	   % model.use_velocitymatrix_file and filecaching mode 2
	   if model.filecache_velocity_matrixfile_extract == 2;
	     model.velocity_matrixfile = ... 
		 cache_velocity_matrixfile_extract(model,PP(1,:),PP(2,:),...
						   'dirichlet_bnd');
	   end;
	   
	   ff = conv_flux(model, repmat(Udir(:),model.flux_quad_degree,1), ...
			  PP(1,:), PP(2,:) );
	   Fdir.Fx = edge_quad_eval_mean(grid,elids,edgeids,...
				      model.flux_quad_degree,ff.Fx);
	   Fdir.Fy = edge_quad_eval_mean(grid,elids,edgeids,...
				      model.flux_quad_degree,ff.Fy);
      end;
      
      %      num_flux.G(dir_NB_ind) = grid.S(dir_NB_ind) .* ...
      %	  (grid.Nx(dir_NB_ind) .* Fdir.Fx + ...
      %	   grid.Ny(dir_NB_ind) .* Fdir.Fy);
      
      % better solution: Lxf with "ghost-cells" 
      num_flux.G(dir_NB_ind) = 1/2 * grid.EL(dir_NB_ind) .* ...
	  (grid.NX(dir_NB_ind).* (flux_mat.Fx(dir_NB_ind) + Fdir.Fx) + ...
	   grid.NY(dir_NB_ind).* (flux_mat.Fy(dir_NB_ind) + Fdir.Fy)) ...
	  + 1/(2 * lambda) * grid.EL(dir_NB_ind) .* ...
	  (UU(dir_NB_ind) - Udir);     
      
    end; % of Dirichlet treatment
    
  %%%%%%%%%%%%%%%%%%%  convective numerical flux: Engquist Osher %%%%%%%%%%%%%%
  elseif isequal(model.name_convective_num_flux,'engquist-osher') ...
	|| isequal(model.name_convective_num_flux,'enquist-osher')
    
    % Engquist-Osher-Flux in 2D:
    %   g_jl(u,v) = |e_jl| * (c^+_jl(u) + c^-_jl(v))
    %   with
    %   c^+_jl(u) = c_jl(0) + int_0^u   max(c_jl'(s),0) ds
    %   c^-_jl(u) =           int_0^u   min(c_jl'(s),0) ds
    %   c_jl(u) = f(u) * n_jl
    %
    %   if we assume, that 
    %         1. f'(u)*n does not change sign on [0,u] and
    %         2. f(0) = 0
    %  then
    %    this simplifies to
    %         c^+_jl(u) = / f(u)*n_jl    if   f'(u)*n_jl > 0
    %                     \ 0            otherwise
    %         c^-_jl(v) = / f(v)*n_jl    if   f'(v)*n_jl < 0
    %                     \ 0            otherwise
    
    
    if model.divclean_mode>0
      error('divergence cleaning not implemented for engquist-osher!')
    end;
    
    % evaluate flux matrix
    % FxNB and FyNB and UNB are computed, where possible:
    %    both are valid in real_NB_ind
    
    % matrix evaluating the flux by quadratures over all edges
    [flux_mat, lambda] = fv_conv_flux_matrix(model,model_data,U);
    
    % matrix with flux values, the neighbour values U(j) inserted
    % values must be reflected across all edges
    % i.e. GNB(i,j) = flux.G(NBI(i,j),INB(i,j))
    % take care of boundary-indices!!
    FxNB = nan * ones([grid.nelements,grid.nneigh]);
    FyNB = nan * ones([grid.nelements,grid.nneigh]);
    %[i,j] = ind2sub(size(flux_mat.G),NB_real_NB_ind);
    [i,j] = ind2sub([grid.nelements, grid.nneigh],NB_real_NB_ind);
    Fx = reshape(flux_mat(1,:,:),size(FxNB));
    Fy = reshape(flux_mat(2,:,:),size(FxNB));
    FxNB(real_NB_ind) = Fx(NB_real_NB_ind);
     FyNB(real_NB_ind) = Fy(NB_real_NB_ind);
     
    
    % evaluate flux matrix derivative (average over edges)
    flux_derivative_mat = fv_conv_flux_linearization_matrix(model,model_data,U);
    % produced fields: Vx, Vy
    
    % matrix with flux values, the neighbour values U(j) inserted
    % values must be reflected across all edges
    % i.e. GNB(i,j) = flux.G(NBI(i,j),INB(i,j))
    % take care of boundary-indices!!
    VxNB = nan * ones(size(flux_derivative_mat.Vx));
    %    [i.j] = ind2sub(size(flux_mat.G),real_NB_ind);
    VxNB(real_NB_ind) = flux_derivative_mat.Vx(NB_real_NB_ind);
    
    VyNB = nan * ones(size(flux_derivative_mat.Vy));
    %    [i.j] = ind2sub(size(flux_mat.G),real_NB_ind);
    VyNB(real_NB_ind) = flux_derivative_mat.Vy(NB_real_NB_ind);
    
    if ~isempty(dir_NB_ind>0)
      % ordinary flux evaluation by quadratures
      %	Fdir = conv_flux(model, Udir, Xdir, Ydir);
      Xdir = grid.ECX(dir_NB_ind);
      Ydir = grid.ECY(dir_NB_ind);
      Udir = model.dirichlet_values_ptr([Xdir,Ydir],model);
      [elids, edgeids] = ind2sub(size(grid.VI),dir_NB_ind);
      PP = edge_quad_points(grid,elids, edgeids,model.flux_quad_degree);

      ff = conv_flux(model, repmat(Udir(:),model.flux_quad_degree,1), ...
		     PP(1,:), PP(2,:) );
      Fdir.Fx = edge_quad_eval_mean(grid,elids,edgeids,...
				    model.flux_quad_degree,ff.Fx);
      Fdir.Fy = edge_quad_eval_mean(grid,elids,edgeids,...
				    model.flux_quad_degree,ff.Fy); 
      
      ff_lin = conv_flux_linearization(model, repmat(Udir(:),...
					      model.flux_quad_degree,...
					      1), ...
				       PP(1,:), PP(2,:) );
      Fdir_derivative.Vx = edge_quad_eval_mean(grid,elids,edgeids,...
						 model.flux_quad_degree,...
						 ff_lin.Vx);
      Fdir_derivative.Vy = edge_quad_eval_mean(grid,elids,edgeids,...
						 model.flux_quad_degree,...
						 ff_lin.Vy); 
      
      FxNB(dir_NB_ind) = Fdir.Fx(:); 
      FyNB(dir_NB_ind) = Fdir.Fy(:); 
      VxNB(dir_NB_ind) = Fdir_derivative.Vx(:); 
      VyNB(dir_NB_ind) = Fdir_derivative.Vy(:); 
      
    end;
    
    % the following nicely maintains nan in Neumann-boundary edges
    c_jl_u_deri = (grid.NX.* flux_derivative_mat.Vx + ...
		   grid.NY.* flux_derivative_mat.Vy);    
    I_c_jl_u_deri_pos = find (c_jl_u_deri > 0);
    I_c_jl_u_deri_nonpos = c_jl_u_deri <= 0;
    
    c_jl_v_deri = (grid.NX.* VxNB + grid.NY.* VyNB);    
    I_c_jl_v_deri_neg = find(c_jl_v_deri < 0);
    I_c_jl_v_deri_nonneg = c_jl_v_deri >= 0;
    
    c_jl_u_plus = nan * ones(size(grid.NBI)); % matrix of Nans
    c_jl_u_plus(I_c_jl_u_deri_pos) = ...
	grid.NX(I_c_jl_u_deri_pos).* Fx(I_c_jl_u_deri_pos) + ...
	grid.NY(I_c_jl_u_deri_pos).* Fy(I_c_jl_u_deri_pos);
    c_jl_u_plus(I_c_jl_u_deri_nonpos) = 0;
    
    c_jl_v_minus = nan * ones(size(grid.NBI));
    c_jl_v_minus(I_c_jl_v_deri_neg) = ...
	grid.NX(I_c_jl_v_deri_neg).* FxNB(I_c_jl_v_deri_neg) + ...
	grid.NY(I_c_jl_v_deri_neg).* FyNB(I_c_jl_v_deri_neg);
    c_jl_v_minus(I_c_jl_v_deri_nonneg) = 0;
 
    % add both components for final flux in all inner edges
    num_flux.G = grid.EL.*(c_jl_u_plus + c_jl_v_minus);
    
    % Lipschitz konstant Lg
    %       |g_ij(w,v) - g_ij(w',v')| <= Lg |S_ij| ( |w'-w| + |v'-v| ) 
    
    num_flux.Lg = 1/ lambda;
    
    %%%%%%%%%%%%%%%%%%%  convective numerical flux: None %%%%%%%%%%%%%%  
  elseif  isequal(model.name_convective_num_flux,'none')
    % do nothing;
    num_flux.Lg = 0;    
  else
    error ('convective numerical flux unknown');
  end;
  
%keyboard;
%| \docupdate 
