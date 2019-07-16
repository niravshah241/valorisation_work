function [L_E_conv,bdir_E_conv] = fv_operators_conv_explicit_engquist_osher(...
                                    model,model_data,U,NU_ind)
%function [L_E_conv,bdir_E_conv] = fv_operators_conv_explicit_engquist_osher(...
%                                   model,model_data[,U,NU_ind])
% computes convection contribution to finite volume time evolution matrices,
% <b> or their Frechet derivative </b>
%
% This function computes a convection operator `L_{\text{eo}}` and a
% corresponding offset vector `b_{\text{eo}}` that can be used by
% fv_operators_implicit_explicit() to build evolution matrices for a finite
% volume time step
% `L_I U^{k+1} = L_E U^k + b_E + b_I`.
%
% @verbatim
% (L_E_conv)il = 1/|Ti| *
%                    (sum j (NB(i) cup NB_dir(i): v_ij*n_ij>=0 )
%                 |S_ij| (v(c_ij)*n_ij))    for l=i
%                1/|Ti| * |S_il| ( - v(c_il)*n_il)
%                                           for l NB(i) and v_il * n_il < 0
%                       0                   else
% @endverbatim
%
% ``
% (L_{E,conv})_{il} =
%     \begin{cases}
%              \frac{1}{|T_i|}
%              \sum_{ j \in \{ NB(i) \cup NB_{dir}(i): v_{ij} \cdot n_{ij}>=0 \} }
%                 |S_{ij}| (v(c_{ij}) \cdot n_{ij})
%                     & \text{for } l=i \\
%               - \frac{1}{|T_i|} |S_{il}|  v(c_{il}) \cdot n_{il}
%                     & \text{for } l \in NB(i)
%                       \text{ and }v_{il} \cdot n_{il} < 0 \\
%               0     &  \text{else}
%     \end{cases}
% ``
%
% The analytical term inspiring this operator looks like
% ` v \cdot \nabla u `
% or in the non-linear case, where the Frechet derivative is computed
% ` \nabla \cdot f(u). `
% Here, `v` is a space dependent velocity field and
% `f` some smooth function
% in `C^2(\mathbb{R}, \mathbb{R}^d)`.
%
% The result are a sparse matrix 'L_E_conv' and an offset vector
% 'bdir_E_conv', the latter containing dirichlet value contributions
%
% See also: fv_num_conv_flux_engquist_osher()
%
% Return values:
%  L_E_conv : sparse matrix `L_{\text{eo}}`
%  bdir_E_conv : offset vector `b_{\text{eo}}`
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

% Bernard Haasdonk 13.7.2006

% determine affine_decomposition_mode as integer
decomp_mode = model.decomp_mode;

grid = [];
if ~isempty(model_data)
  grid = model_data.grid;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% engquist osher %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

  % ATTENTION: by using u=1, we get the velocity field
  % Of course, this only works, when the conv_flux function is linear
  if model.flux_linear
    flux_mat = fv_conv_flux_matrix(model,model_data,ones(n,1));
  else
    flux_mat_sign = fv_conv_flux_matrix(model, model_data, ones(n,1));
    % non-linear case: use f' for flux matrix generation
    flux_mat = fv_conv_flux_matrix(model, model_data, U, model.conv_flux_derivative_ptr);
  end
  %  disp('halt after flux_mat')
  %  keyboard;
else % decomp_mode == 2
  flux_mat = fv_conv_flux_matrix(model,[],[]);
end;

if model.debug && decomp_mode == 0
  disp('test affine decomposition of flux_matrix:');
  test_affine_decomp(@fv_conv_flux_matrix,1,1,model,model_data, ...
                     ones(n,1));
  disp('OK');
end;

if decomp_mode == 2
  L_E_conv = flux_mat(:);
elseif decomp_mode == 0
  V.Vx = reshape(flux_mat(1,:,:),size(grid.ECX));
  V.Vy = reshape(flux_mat(2,:,:),size(grid.ECX));

  %%%%%%%% explicit matrix:
  % has entries
  % (L_E_conv)il = 1/|Ti| *
  %                    (sum j (NB(i) cup NB_dir(i): v_ij*n_ij>=0 )
  %                 |S_ij| (v(c_ij)*n_ij))         for l=i
  %                1/|Ti| * |S_il| ( - v(c_il)*n_il)
  %                        for l NB(i) and v_il * n_il < 0
  %                                                       0       else

  % compute products
  vn = zeros(size(grid.ESX));
  vn(real_or_dir_NB_ind) = grid.EL(real_or_dir_NB_ind) .* ...
      (V.Vx(real_or_dir_NB_ind).*grid.NX(real_or_dir_NB_ind)+ ...
       V.Vy(real_or_dir_NB_ind).*grid.NY(real_or_dir_NB_ind));
  %    la = zeros(size(grid.ESX));
  %    la(real_or_dir_NB_ind) =         0.5 * grid.EL(real_or_dir_NB_ind) * ...
  %     ( 1 / lambda);
  if model.flux_linear
    vn_sign = vn;
  else
    vn_sign = zeros(size(grid.ESX));
    V_sign.Vx = reshape(flux_mat(1,:,:),size(grid.ECX));
    V_sign.Vy = reshape(flux_mat(2,:,:),size(grid.ECX));
    vn_sign(real_or_dir_NB_ind) = grid.EL(real_or_dir_NB_ind) .* ...
      (V_sign.Vx(real_or_dir_NB_ind).*grid.NX(real_or_dir_NB_ind)+ ...
       V_sign.Vy(real_or_dir_NB_ind).*grid.NY(real_or_dir_NB_ind));
  end

  % diagonal entries
  fi_pos = find(vn_sign>0);
  fi_neg = find(vn_sign<0);
  vn_pos = zeros(size(grid.ESX));
  vn_neg = zeros(size(grid.ESX));
  vn_pos(fi_pos) = vn(fi_pos);
  vn_neg(fi_neg) = vn(fi_neg);
  L_E_diag = sparse(1:n,1:n, grid.Ainv(:) .* sum(vn_pos,2));

  % off-diagonal entries (nonpositive):
  [i,dummy ]= ind2sub(size(grid.ESX),real_NB_ind);
  L_E_offdiag = sparse(i,grid.NBI(real_NB_ind),grid.Ainv(i) .* ...
                       vn_neg(real_NB_ind), n,n);
  L_E_conv = L_E_diag + L_E_offdiag;

  if ~isempty(NU_ind)
    L_E_conv = L_E_conv(NU_ind, :);
  end
  bdir_E_conv = 0;
  return;

elseif decomp_mode==1

  % components as
  % v-decomposition

  Q_v = length(flux_mat);
  L_E_conv = cell(Q_v,1);
  % auxiliary quantities required some times
  [i,dummy ]= ind2sub(size(grid.ESX),real_NB_ind);

  % velocity components
  for q = 1:Q_v;
    V.Vx = reshape(flux_mat{q}(1,:,:),size(grid.ECX));
    V.Vy = reshape(flux_mat{q}(2,:,:),size(grid.ECX));
    % compute products
    vn = zeros(size(grid.ESX));
    vn(real_or_dir_NB_ind) = grid.EL(real_or_dir_NB_ind) .* ...
        (V.Vx(real_or_dir_NB_ind).*grid.NX(real_or_dir_NB_ind)+ ...
         V.Vy(real_or_dir_NB_ind).*grid.NY(real_or_dir_NB_ind));

    % diagonal entries
    fi_pos = find(vn>0);
    fi_neg = find(vn<0);
    vn_pos = zeros(size(grid.ESX));
    vn_neg = zeros(size(grid.ESX));
    vn_pos(fi_pos) = vn(fi_pos);
    vn_neg(fi_neg) = vn(fi_neg);

    % diagonal entries
    L_E_diag = sparse(1:n,1:n,   grid.Ainv(:) .* sum(vn_pos,2));
    % off-diagonal entries (nonpositive):
    L_E_offdiag = sparse(i,grid.NBI(real_NB_ind),grid.Ainv(i) .* ...
                         (vn_neg(real_NB_ind)), n,n);
    L_E_conv{q} = L_E_diag + L_E_offdiag;

  end;
end;

%%%%%%%% dirichlet-offset-vector:

% the following is a very rough discretization, should be replaced
% by higher order quadratures sometime, as the flux-matrix also is
% evaluated with higher order
%
% (bdir_E_conv)_i = - 1/|T_i| * ...
%    sum_j Ndir_vn_negative(i) ( |S_ij| v(c_ij)*n_ij u_dir(c_ij,t)
if decomp_mode == 2 % decomp_mode == 2 -> coefficients
  Q_v = length(flux_mat);
  Udir = model.dirichlet_values_ptr([],model);
  Q_Udir = length(Udir);
  bdir_E_conv = zeros(Q_Udir * Q_v,1);
  for q1 = 1:Q_Udir
    for q2 = 1:Q_v
      bdir_E_conv((q1-1)*(Q_v)+ q2) = Udir(q1)*flux_mat(q2);
    end;
  end;


elseif decomp_mode == 0
  if ~isempty(dir_NB_ind > 0)
    % evaluate dirichlet values at edge midpoints at t
    Xdir = grid.ECX(dir_NB_ind);
    Ydir = grid.ECY(dir_NB_ind);
    Udir = model.dirichlet_values_ptr([Xdir(:),Ydir(:)],model);

    if model.debug
      disp('test affine decomposition of dirichlet values:');
      test_affine_decomp(model.dirichlet_values_ptr,...
                         1,2,[Xdir(:), Ydir(:)],model);
      disp('OK');
    end;

    val2 = zeros(size(grid.ESX));
%    vn_dir_nb = zeros(length(dir_NB_ind),1);
    vn_dir_nb = ...
        V.Vx(dir_NB_ind).*grid.NX(dir_NB_ind) + ...
        V.Vy(dir_NB_ind).*grid.NY(dir_NB_ind);
    fi_neg = find(vn_dir_nb<0);
    val2(dir_NB_ind(fi_neg))= grid.EL(dir_NB_ind(fi_neg)).* ...
        (vn_dir_nb(fi_neg)) .* Udir(fi_neg);
    bdir_E_conv = - grid.Ainv(:).* sum(val2,2);
%    save('ws_complete');
  else
    bdir_E_conv = zeros(grid.nelements,1);
  end;

else % decomp_mode == 1
  if ~isempty(dir_NB_ind > 0)
    % for each udir-component compute Q_v components

    % evaluate dirichlet values at edge midpoints at t
    Xdir = grid.ECX(dir_NB_ind);
    Ydir = grid.ECY(dir_NB_ind);
    Udir = model.dirichlet_values_ptr([Xdir(:),Ydir(:)],model);

    Q_Udir = length(Udir);
    bdir_E_conv = cell(Q_Udir * Q_v,1);

%    vn_dir_nb = zeros(length(dir_NB_ind),1);
    for q1 = 1:Q_Udir
      for q2 = 1:Q_v
        val2 = zeros(size(grid.ESX));
        Fx = reshape(flux_mat{q2}(1,:,:),size(grid.ECX));
        Fy = reshape(flux_mat{q2}(2,:,:),size(grid.ECX));
        vn_dir_nb = ...
            Fx(dir_NB_ind).*grid.NX(dir_NB_ind) + ...
            Fy(dir_NB_ind).*grid.NY(dir_NB_ind);
        fi_neg = find(vn_dir_nb<0);
        val2(dir_NB_ind(fi_neg))= grid.EL(dir_NB_ind(fi_neg)).* ...
            (vn_dir_nb(fi_neg)) .* (Udir{q1}(fi_neg));
%       disp(['val2 nonzeros:',num2str(length(find(val2))),...
%            ' q1 = ',num2str(q1),' q2=',num2str(q2)]);
        %         val2(dir_NB_ind)= grid.EL(dir_NB_ind).* ...
        %             (flux_mat{q2}.Fx(dir_NB_ind).*grid.NX(dir_NB_ind) + ...
        %              flux_mat{q2}.Fy(dir_NB_ind).*grid.NY(dir_NB_ind)) .*Udir{q1};
        bdir_E_conv{(q1-1)*Q_v+q2} = - grid.Ainv(:).* sum(val2,2);
      end;
    end;
%    save('ws_components');
  else
    % still produce dummy components as the availability of dir
    % boundary cannot be detected in online phase
    tmodel = model;
    tmodel.decomp_mode = 2;
    Udir = model.dirichlet_values_ptr([],tmodel);
    % for each combination of Udir and diffusivity component
    % perform identical computation as above
    Q_Udir = length(Udir);
    bdir_E_conv = cell(Q_Udir * Q_v,1);
  end;
end;
