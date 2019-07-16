function num_flux = fv_num_conv_flux_engquist_osher(model,model_data,U,NU_ind)
%function num_flux = fv_num_conv_flux_engquist_osher(model,model_data,U,NU_ind)
% Function computing a numerical convective Engquist-Osher
% flux matrix.
%
% Engquist-Osher-Flux in 2D:
%   ``g_{jl}(u,v) = |e_{jl}| \cdot (c^+_{jl}(u) + c^-_{jl}(v))``
%   with
% ``c^+_{jl}(u) = c_{jl}(0) + \int_0^u   \max(c_{jl}'(s),0) ds``
% ``c^-_{jl}(u) =             \int_0^u   \min(c_{jl}'(s),0) ds``
% ``c_{jl}(u)   = f(u) \cdot n_{jl}``
%
%   if we assume, that
%         -# `f'(u)\cdot n` does not change sign on `[0,u]` and
%         -# `f(0) = 0`
%  then
%    this simplifies to
%     ``c^+_{jl}(u) = \left\{\begin{array}{rl}
%         f(u)\cdot n_{jl} & \mbox{if }f'(u)\cdot n_{jl} > 0 \\
%         0                & \mbox{otherwise}
%       \end{array}\right.``
%     ``c^-_{jl}(v) = \left\{\begin{array}{rl}
%         f(v)\cdot n_{jl} & \mbox{if }f'(v)\cdot n_{jl} < 0 \\
%         0                & \mbox{otherwise}.
%       \end{array}\right.``
%
% Dirichlet-boundary treatment is performed, but
% Neumann-values are set to 'nan'. They must be performed by
% the calling function.
%
% fields of num_flux:
%   Lg: lipschitz constant satisfying
%       ``|g_{ij}(w,v) - g_{ij}(w',v')| \leq Lg |S_{ij}| ( |w'-w| + |v'-v| ).``
%   G:  the matrix of numerical flux values across the edges.
%       boundary treatment is performed as follows:
%         - in dirichlet-boundary points, the
%         neighbour-values are simply set to this value and
%         the flux is computed.
%         - Neumann-treatment is not performed here
%         .
%
% We assume for the convective flux function:
%         -# `f'(u)\cdot n` does not change sign on interval `(u,0)`
%         -# `f(0) = 0`
%
% Additionally, a velocity function must be available in order to
% compute the derivative of the convection term.
%

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
  
  if ~isempty(model.flux_quad_degree) 
    model.flux_quad_degree = 1;
  end;
  
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
  % UNB = nan * ones(size(grid.NBI));
  % UNB(real_NB_ind) = U(grid.NBI(real_NB_ind));
  UU = repmat(U,1,size(grid.NBI,2));
    
  %%%%%%%%%%%%%%%%%%%  convective numerical flux: Engquist Osher %%%%%%%%%%%%%%

  if model.flux_linear && model.divclean_mode>0
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
    PP = edge_quad_points(grid,elids,edgeids,model.flux_quad_degree);

    ff = model.conv_flux_ptr(PP,...
           repmat(Udir(:),model.flux_quad_degree,1), ...
           model );
    Fdir.Fx = edge_quad_eval_mean(grid,elids,edgeids,...
      model.flux_quad_degree,ff(:,1));
    Fdir.Fy = edge_quad_eval_mean(grid,elids,edgeids,...
      model.flux_quad_degree,ff(:,2)); 

    ff_lin = model.conv_flux_derivative_ptr(PP,...
               repmat(Udir(:),model.flux_quad_degree, 1), ...
               model);
    Fdir_derivative.Vx = edge_quad_eval_mean(grid,elids,edgeids,...
      model.flux_quad_degree,...
      ff_lin(:,1));
    Fdir_derivative.Vy = edge_quad_eval_mean(grid,elids,edgeids,...
      model.flux_quad_degree,...
      ff_lin(:,2)); 

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

%  if grid.nelements == 200
%    disp('halt in fv_num_conv_flux_engquist_osher');
%    keyboard;
%  end;

%keyboard;
