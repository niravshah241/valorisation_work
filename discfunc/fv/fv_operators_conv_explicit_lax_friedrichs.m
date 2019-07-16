function [L_E_conv,bdir_E_conv] = fv_operators_conv_explicit_lax_friedrichs(...
    model,model_data,U,NU_ind)
%function [L_E_conv,bdir_E_conv] = fv_operators_conv_explicit_lax_friedrichs(...
%   model,model_data[,U,NU_ind])
% computes convection contribution to finite volume time evolution matrices,
% <b> or their Frechet derivative </b>
%
% This function computes a convection operator `L_{\text{lf}}` and a
% corresponding offset vector `b_{\text{lf}}` that can be used by
% fv_operators_implicit_explicit() to build evolution matrices for a finite
% volume time step
% `L_I U^{k+1} = L_E U^k + b_E + b_I`.
%
% The analytical term inspiring this operator looks like
% `` v \cdot \nabla u ``
% or in the non-linear case, where the Frechet derivative is computed
% `` \nabla f(u). ``
% Here, `v` is a space dependent velocity field and `f` some smooth function
% in `C^2(\mathbb{R}^d, \mathbb{R})`. The fluxes can be controlled with the
% 'model.conv_flux_ptr' fields. The first case is default usage of this
% operator with the field 'params.conv_flux_ptr' pointing to
% conv_flux_linear() where the velocity field `v` is given through
% 'velocity_ptr'.
%
% The Lax-Friedrichs flux is used with fluxes:
%   ``g_{jl}(u,v) = \frac12 (f(u,x_{jl}) + f(v,x_{jl})) \cdot \nu_{jl} -
%             \frac{|S_{jl}|}{2 \lambda_{jl}} (v-u)``
%   with `\lambda_{jl} := \frac{\Delta t}{|m_j - m_l|}`.
%
% See also: fv_num_conv_flux_lax_friedrichs().
%
% Required fields of model:
%  lxf_lambda : scalar `\lambda = \sup \lambda_{jl}` controlling the
%               artificial diffusion added to the flux.
%  dirichlet_values_ptr : function pointer for dirichlet function
%                         `u_{\text{dir}}`
%  flux_linear : flag indicating wether the flux function `f` is linear. If
%                this flag is set to 'false' and the decomp mode is set to
%                'complete' the Frechet derivative `DL_{\text{lf}}|_{u}` is
%                returned.
%  conv_flux_ptr : function pointer specifying the function `f`.
%
% Optional fields of model:
%  conv_flux_derivative_ptr : function pointer specifying the derivative of
%                             the convection function `f`. This is only needed
%                             in case we want to compute the Frechet
%                             derivative.
%
% Return values:
%  L_E_conv : sparse matrix `L_{\text{lf}}` or `DL_{\text{lf}}|_{u}` if flux
%             is non-linear.
%  bdir_E_conv : offset vector `b_{\text{lf}}`
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

% determine affine_decomposition_mode as integer
decomp_mode = model.decomp_mode;

grid = [];
if ~isempty(model_data)
  grid = model_data.grid;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% lax-friedrichs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if decomp_mode < 2
  n = grid.nelements;  
  % determine index sets of boundary types
  real_NB_ind = find(grid.NBI>0);
  %NBIind = grid.NBI(real_NB_ind);
  %INBind = grid.INB(real_NB_ind);
  dir_NB_ind =find(grid.NBI == -1);
  real_or_dir_NB_ind = [real_NB_ind(:); dir_NB_ind(:)];
  if model.verbose >= 29
    disp(['found ',num2str(length(real_NB_ind)),' non-boundary edges.'])
    disp(['found ',num2str(length(dir_NB_ind)),' dirichlet-boundary edges.'])
  end;

  if nargin < 4
    NU_ind = [];
  end

  lambda   = model.lxf_lambda;
  if model.flux_linear
    % by using u=1, we get the velocity field
    flux_mat = fv_conv_flux_matrix(model,model_data,ones(n,1));
  else
    % non-linear case: use f' for flux matrix generation
    old_conv_flux_ptr = model.conv_flux_ptr;
    model.conv_flux_ptr = model.conv_flux_derivative_ptr;
    flux_mat = fv_conv_flux_matrix(model, model_data, U);
    model.conv_flux_ptr = old_conv_flux_ptr;
  end
  %  disp('halt after flux_mat')
  %  keyboard;
else % decomp_mode == 2
  flux_mat = fv_conv_flux_matrix(model,[],[]);
end;

if decomp_mode == 0
  V.Vx = reshape(flux_mat(1,:,:),size(grid.ECX));
  V.Vy = reshape(flux_mat(2,:,:),size(grid.ECX));
  
  %%%%%%%% explicit matrix:
  % has entries 
  % (L_E_conv)il = 1/|Ti| *
  %                    (sum j (NB(i) cup NB_dir(i)) 
  %                 1/2 |S_ij| (v(c_ij)*n_ij + 1/lambda))         for l=i   
  %              1/|Ti| * 1/2 |S_il|  (v(c_il)*n_il - 1/lambda)   for l NB(i)
  %                                                       0       else 
  
  % compute products    
  vn = zeros(size(grid.ESX));
  vn(real_or_dir_NB_ind) = 0.5 * grid.EL(real_or_dir_NB_ind) .* ...
    (V.Vx(real_or_dir_NB_ind).*grid.NX(real_or_dir_NB_ind)+ ...
    V.Vy(real_or_dir_NB_ind).*grid.NY(real_or_dir_NB_ind));
  la = zeros(size(grid.ESX));
  la(real_or_dir_NB_ind) = 	0.5 * grid.EL(real_or_dir_NB_ind) * ...
    ( 1 / lambda);
  
  % diagonal entries
  L_E_diag = sparse(1:n,1:n,   grid.Ainv(:) .* sum(vn+la,2));
  
  % off-diagonal entries (nonpositive):
  [i,dummy ]= ind2sub(size(grid.ESX),real_NB_ind);
  L_E_offdiag = sparse(i,grid.NBI(real_NB_ind),grid.Ainv(i) .* ...
		 (vn(real_NB_ind) - la(real_NB_ind)), n,n);
  L_E_conv = L_E_diag + L_E_offdiag;

  % check for stability: off-diagonal entries must be nonpositive
  % (mult with -dt gives the coefficients Theta_ij in the convex combination
  %  u_i^(k+1) = (1-sum Theta_ij)u_i^k + sum Theta_ij u_j^k)
  [i,j] = find(L_E_offdiag>0);
  if model.verbose>=9
    if ~isempty(i)
      error('error: lambda chosen too large! non stable LxF-scheme!');
    end;
  end;
  if ~isempty(NU_ind)
    L_E_conv = L_E_conv(NU_ind, :);
  end

  %    disp('debug halt:')
  %    keyboard;

  % check row-sum of L_E_conv: must be positive! 
  %   -> no, only for non-neuman-elements
  %if model.verbose>=10
  %  if ~isempty(find(sum(L_E_conv,2)<0))
  %	disp('row sum of explicit matrix is smaller than zero!!')
  %	% derivation yields:
  %	% sum(L_E_conv,2)_i = T1(i) + T2(i) 
  %	% with
  %	% T1(i) = 1/|Ti| sum_ (j NB(i) or NB_dir(i)) |Sij|(v_ij n_ij)
  %	% T2(i) = 1/|Ti| sum_ (j NB_dir(i)) 0.5* |Sij|(1/lambda - v_ij n_ij)
  %	m1 = zeros(size(grid.Sx));
  %	m2 = zeros(size(grid.Sx));
  %	m1(real_or_dir_NB_ind) = grid.S(real_or_dir_NB_ind).* ...
  %	    (V.Vx(real_or_dir_NB_ind).*grid.Nx(real_or_dir_NB_ind)+ ...
  %	     V.Vy(real_or_dir_NB_ind).*grid.Ny(real_or_dir_NB_ind));
  %	m2(dir_NB_ind) = 0.5 * grid.S(dir_NB_ind) .* ...
  %	    (1/lambda - ...
  %	     V.Vx(dir_NB_ind).*grid.Nx(dir_NB_ind)- ...
  %	     V.Vy(dir_NB_ind).*grid.Ny(dir_NB_ind));
  %	T1 = sum(m1,2).*grid.Ainv(:);
  %	T2 = sum(m2,2).*grid.Ainv(:);
  %	disp('check T1+T2 = sum(L_E_conv,2) and nonnegativity of T1, T2');
  %	keyboard;
  %      end;
  %    end;
  
elseif decomp_mode==1
  if ~model.flux_linear
    error('flux needs to be linear')
  end
  % first component is lambda-term, remaining terms as
  % v-decomposition
  
  Q_v = length(flux_mat);
  L_E_conv = cell(1+Q_v,1);
  % auxiliary quantities required some times
  [i,dummy ]= ind2sub(size(grid.ESX),real_NB_ind);
  
  % lambda-component    
  la = zeros(size(grid.ESX));
  la(real_or_dir_NB_ind) = 	0.5 * grid.EL(real_or_dir_NB_ind) * ...
    ( 1 / lambda);    
  % diagonal entries
  L_E_diag = sparse(1:n,1:n,   grid.Ainv(:) .* sum(la,2));    
  % off-diagonal entries (nonpositive):
  L_E_offdiag = sparse(i,grid.NBI(real_NB_ind),grid.Ainv(i) .* ...
		 ( - la(real_NB_ind)), n,n);
  L_E_conv{1} = L_E_diag + L_E_offdiag;
  
  % velocity components
  for q = 1:Q_v;
    V.Vx = reshape(flux_mat{q}(1,:,:), size(grid.ECX));
    V.Vy = reshape(flux_mat{q}(2,:,:), size(grid.ECX));
    % compute products    
    vn = zeros(size(grid.ESX));
    vn(real_or_dir_NB_ind) = 0.5 * grid.EL(real_or_dir_NB_ind) .* ...
      (V.Vx(real_or_dir_NB_ind).*grid.NX(real_or_dir_NB_ind)+ ...
      V.Vy(real_or_dir_NB_ind).*grid.NY(real_or_dir_NB_ind));
    % diagonal entries
    L_E_diag = sparse(1:n,1:n,   grid.Ainv(:) .* sum(vn,2));    
    % off-diagonal entries (nonpositive):
    L_E_offdiag = sparse(i,grid.NBI(real_NB_ind),grid.Ainv(i) .* ...
		   (vn(real_NB_ind)), n,n);
    L_E_conv{q+1} = L_E_diag + L_E_offdiag;
  end;
else % decomp_mode == 2 -> coefficients
  L_E_conv = [1; flux_mat(:)];
end;

%%%%%%%% dirichlet-offset-vector:
% (bdir_E_conv)_i = - 1/|T_i| * ...
%    sum_j Ndir(i) (1/2 |S_ij| (v(c_ij)*n_ij-1/lambda) u_dir(c_ij,t))
if decomp_mode == 0
  if ~isempty(dir_NB_ind > 0)
    % evaluate dirichlet values at edge midpoints at t
    Xdir = grid.ECX(dir_NB_ind);
    Ydir = grid.ECY(dir_NB_ind);
    Udir = model.dirichlet_values_ptr([Xdir,Ydir],model);
    
    val2 = zeros(size(grid.ESX));
    val2(dir_NB_ind)= 0.5 * grid.EL(dir_NB_ind).* ... 
      (V.Vx(dir_NB_ind).*grid.NX(dir_NB_ind) + ...
      V.Vy(dir_NB_ind).*grid.NY(dir_NB_ind)-1/lambda) .*Udir; 
    bdir_E_conv = - grid.Ainv(:).* sum(val2,2);      
  else
    bdir_E_conv = zeros(grid.nelements,1);
  end;
elseif decomp_mode == 1
  if ~isempty(dir_NB_ind > 0)
    % for each u-component compute 1+Q_v components
    
    % evaluate dirichlet values at edge midpoints at t
    Xdir = grid.ECX(dir_NB_ind);
    Ydir = grid.ECY(dir_NB_ind);
    Udir = model.dirichlet_values_ptr([Xdir, Ydir], model);
    
    Q_Udir = length(Udir);
    bdir_E_conv = cell(Q_Udir * (Q_v + 1),1);
    
    val2 = zeros(size(grid.ESX));
    for q1 = 1:Q_Udir
      % first component is lambda
      val2(dir_NB_ind)= 0.5 * grid.EL(dir_NB_ind) * ... 
        (-1/lambda) .*Udir{q1}; 
      bdir_E_conv{(q1-1)*(Q_v +1)+1} = - grid.Ainv(:).* sum(val2,2);      
      % remaining are velocity-components
      for q2 = 1:Q_v
        F.Fx = reshape(flux_mat{q2}(1,:,:),size(grid.ECX));
        F.Fy = reshape(flux_mat{q2}(2,:,:),size(grid.ECX));
        val2(dir_NB_ind)= 0.5 * grid.EL(dir_NB_ind).* ... 
          (F.Fx(dir_NB_ind).*grid.NX(dir_NB_ind) + ...
           F.Fy(dir_NB_ind).*grid.NY(dir_NB_ind)) .*Udir{q1}; 
        bdir_E_conv{(q1-1)*(Q_v+1)+1+q2} = - grid.Ainv(:).* sum(val2,2);
      end;
    end;
  else
    % still produce dummy components as the availability of dir
    % boundary cannot be detected in online phase 
    tmodel = model;
    tmodel.decomp_mode = 2;
    Udir = model.dirichlet_values_ptr([],tmodel);
    % for each combination of Udir and diffusivity component
    % perform identical computation as above
    Q_Udir = length(Udir);
    bdir_E_conv = cell(Q_Udir * (Q_v+1),1);
  end; 
  
else % decomp_mode == 2 -> coefficients
  Q_v = length(flux_mat);
  Udir = model.dirichlet_values_ptr([],model);
  Q_Udir = length(Udir);
  bdir_E_conv = zeros(Q_Udir * (Q_v + 1),1);      
  for q1 = 1:Q_Udir
    bdir_E_conv((q1-1)*(Q_v +1)+1) = Udir(q1);      
    for q2 = 1:Q_v
      bdir_E_conv((q1-1)*(Q_v+1)+1+q2) = Udir(q1)*flux_mat(q2);      
    end;
  end;
end;


