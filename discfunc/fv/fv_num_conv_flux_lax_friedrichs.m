function num_flux = fv_num_conv_flux_lax_friedrichs(model,model_data,U,NU_ind)
% function num_flux = fv_num_conv_flux_lax_friedrichs(model,model_data,U,NU_ind)
% Function computing a numerical convective Lax-Friedrichs
% flux matrix.
%
% Dirichlet-boundary treatment is performed, but Neumann-values are set to
% 'nan'. They must be performed by the calling function.
%
% Lax-Friedrichs-Flux:
%   `g_{jl}(u,v) = \frac12 (f(u,x_{jl}) + f(v,x_{jl})) \cdot \nu_{jl} -
%             \frac{|S_{jl}|}{2 \lambda_{jl}} (v-u)`
%   or equivalently
%   `g_{jl}(u,v) = \frac12 |S_{jl}| \left( (f(u)+f(v)) \cdot n_{jl} -
%                           \frac{1}{\lambda_{jl}} (v-u) \right)`
%   with `\lambda_{jl} := \frac{\Delta t}{|m_j - m_l|}`.
%
% generated fields of num_flux:
%   Lg: Lipschitz constant satisfying
%       `|g_{ij}(w,v) - g_{ij}(w',v')| \leq Lg |S_{ij}| ( |w'-w|
%       + |v'-v| )`
%   G:  the matrix of numerical flux values across the edges.
%       boundary treatment is performed as follows:
%         - in dirichlet-boundary points, the
%         neighbour-values are simply set to this value and
%         the flux is computed.  
%         - Neumann-treatment is not performed here
%         .
%
% required fields of model:
%   lxf_lambda: numerical diffusion in Lax-Friedrichs flux

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
% NB_real_ind comprises some kind of inverse cell-neighbour relation that
% leads back to the original cell from where we started:
%
% If
%  [i,j]=ind2sub(size(grid.NBI),NB_real_NB_ind(sub2ind(size(grid.NBI),k,l)))
% then
%  grid.NBI(i,j) == k
NB_real_NB_ind = sub2ind(size(grid.NBI),NBIind,INBind);

dir_NB_ind = find(grid.NBI == -1);

% matrix with the neighbouring function values
% UNB = nan * ones(size(flux_mat.Fx));
UNB = nan * ones(size(grid.NBI));
UNB(real_NB_ind) = U(grid.NBI(real_NB_ind));
UU = repmat(U,1,size(grid.NBI,2));

% 'FxNB' and 'FyNB' and 'UNB' are computed, where possible:
%    both are valid in 'real_NB_ind'

% matrix evaluating the flux by quadratures over all edges
[flux_mat,lambda] = fv_conv_flux_matrix(model,model_data,U);

% in case of divergence cleaning, the velocity field must be
% available!! This is obtianed by setting u to 1
if (model.divclean_mode>0)
  Udummy = ones(size(U));
  [flux_mat_Uone,lambda_one] = fv_conv_flux_matrix(model,model_data,Udummy);
end;

% matrix with flux values, the neighbour values U(j) inserted
% values must be reflected across all edges
% i.e. GNB(i,j) = flux.G(NBI(i,j),INB(i,j))
% take care of boundary-indices!!
FxNB = nan * ones([grid.nelements,grid.nneigh]);
FyNB = nan * ones([grid.nelements,grid.nneigh]);
% [i,j] = ind2sub([grid.nelements, grid.nneigh],NB_real_NB_ind);
Fx = reshape(flux_mat(1,:,:),size(FxNB));
Fy = reshape(flux_mat(2,:,:),size(FxNB));
FxNB(real_NB_ind) = Fx(NB_real_NB_ind);
FyNB(real_NB_ind) = Fy(NB_real_NB_ind);

% L1 is the supremum in space of the derivative of the first component
% L2 is the supremum in space of the derivative of the second component
L1 = max(max(abs(Fx)));
L2 = max(max(abs(Fy)));

% Lipschitz-constant of LxF Flux
num_flux.Lg = 1/(2*lambda) + 1/2 * norm([L1, L2]);
% an alternative formula without motivation:
%    num_flux.Lg =  max(max(grid.S)) * model.c * 0.25;

% in the linear case the following worked:
% numflux_mat.G = 1/2 * grid.S .* ...
%               (grid.Nx.* (flux_mat.Vx.*(OU + OUext(grid.NBI))) + ...
%             grid.Ny.* (flux_mat.Vy.*(OU + OUext(grid.NBI)))) ...
%       + 1/(2 * flux_mat.lambda) * grid.S .* (OU - OUext(grid.NBI));

if ~isempty(model.lxf_lambda)
  lambda = model.lxf_lambda;
end;

% in the general nonlinear case:
num_flux.G = 1/2 * grid.EL .* ...
    (grid.NX.* (Fx + FxNB) + ...
     grid.NY.* (Fy + FyNB)) ...
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
  if model.divclean_mode > 0
    if model.verbose>9
      disp('skipping boundary flux eval. due to divergence cleaning!!');
    end;
    % the following is unstable, if U reaches 0!!!!
    %   Fdir.Fx = flux_mat.Fx(dir_NB_ind)./UU(dir_NB_ind).*Udir;
    %   Fdir.Fy = flux_mat.Fy(dir_NB_ind)./UU(dir_NB_ind).*Udir;
    % instead again evaluation of the whole matrix.
    FOnex = reshape(flux_mat_UOne(1,:,:),size(FxNB));
    FOney = reshape(flux_mat_Uone(2,:,:),size(FxNB));
    Fdir = [FOnex(dir_NB_ind),FOney(dir_NB_ind)].*Udir;

  else % ordinary flux evaluation by quadratures
       %        Fdir = conv_flux(model,Udir, Xdir, Ydir);
       [elids, edgeids] = ind2sub(size(grid.VI),dir_NB_ind);
       PP = edge_quad_points(grid,elids, edgeids, ...
                             model.flux_quad_degree);

       % the following is only relevant in case of use of a
       % model.use_velocitymatrix_file and filecaching mode 2
       if isfield(model, 'velocity_matrixfile_extract') && ...
           model.filecache_velocity_matrixfile_extract == 2;
         model.velocity_matrixfile = ...
           cache_velocity_matrixfile_extract(model,PP(:,1),PP(:,2),...
                                           'dirichlet_bnd');
       end;

       ff = model.conv_flux_ptr(PP, repmat(Udir(:),model.flux_quad_degree,1), ...
                                model );
       Fdir.Fx = edge_quad_eval_mean(grid,elids,edgeids,...
                              model.flux_quad_degree,ff(:,1));
       Fdir.Fy = edge_quad_eval_mean(grid,elids,edgeids,...
                              model.flux_quad_degree,ff(:,2));
  end;

  %      num_flux.G(dir_NB_ind) = grid.S(dir_NB_ind) .* ...
  %       (grid.Nx(dir_NB_ind) .* Fdir.Fx + ...
  %        grid.Ny(dir_NB_ind) .* Fdir.Fy);

  % better solution: Lxf with "ghost-cells"
  num_flux.G(dir_NB_ind) = 1/2 * grid.EL(dir_NB_ind) .* ...
      (grid.NX(dir_NB_ind).* (Fx(dir_NB_ind) + Fdir.Fx) + ...
       grid.NY(dir_NB_ind).* (Fy(dir_NB_ind) + Fdir.Fy)) ...
      + 1/(2 * lambda) * grid.EL(dir_NB_ind) .* ...
      (UU(dir_NB_ind) - Udir);

end

