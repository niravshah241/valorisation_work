function [L_I_diff, bdir_I_diff] = fv_operators_diff_implicit_gradient(...
    model, model_data, U, NU_ind)
%function [L_I_diff, bdir_I_diff] = ...
%    fv_operators_diff_implicit_gradient(model,model_data[, U, NU_ind])
% computes diffusion contributions to finite volume time evolution matrices,
% <b> or their Frechet derivative </b>
%
% This function computes a diffusion operator `L_{\text{diff}}` and a
% corresponding offset vector `b_{\text{diff}}` that can be used by
% fv_operators_implicit_explicit() to build evolution matrices for a finite
% volume time step
% `L_I U^{k+1} = L_E U^k + b_E + b_I`.
%
% The analytical term inspiring this operator looks like
% `` - \nabla \cdot \left(d(u) \nabla l(u) \right). ``
% Here, the functions `d(u)` and `l(u)` differ in the numerical
% implementation.
%
% The latter is evaluated in cell centers during flux computation, whereas
% the diffusion functional `d(u)` is evaluated on the cell edges by averaging
% the solutions value there.
%
% See also: fv_num_diff_flux_gradient()
%
% The implemented flux then can be written for inner boundaries as
%   `` g_{ij}(u,v) = d(\{u,v\}) \frac{l(v) - l(u)}{|d_{ij}|} |e_{ij}| ``
% and for Dirichlet boundaries as
%   `` g_{dir,ij}(u) = d(u_{dir}(x_{ij}))
%          \frac{l(u_{dir}(x_{ij})) - u}{|d_{ij}|} |e_{ij}|, ``
% where `\{u,v\}` refers to the average of `u` and `v`.
%
% \note that this implementation only works for trivial functions `l = id`,
% because otherwise the operator would became non-linear in the space
% variable.
%
% If, however, `l(u)` is non-trivial, the flux arguments are simply
% <b>scaled</b> by `l'(u)` giving rise to the Frechet derivative
% `D L_I|_{U}` that can be used in a Newton scheme. The offset vector
% 'bdir_I_diff' can then be ignored, of course!
%
% required fields of model:
%   diffusivity_ptr          : The diffusivity function `d(u)`
%   laplacian_derivative_ptr : The derivative `l'(u)`. See above for
%                              explanation.
%   dirichlet_values_ptr     : The dirichlet function `u_{\text{dir}}`
%
% Return values:
%  L_I_diff       : a sparse matrix with diffusion contributions to `L_I`.
%  bdir_I_diff    : and offset vector containing the Dirichlet value
%                   contributions of the diffusion parts.
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

if nargin < 4
  NU_ind = [];
end

if decomp_mode < 2
  n = grid.nelements;
  % determine index sets of boundary types
  real_NB_ind = find(grid.NBI>0);
  %NBIind = grid.NBI(real_NB_ind);
  %INBind = grid.INB(real_NB_ind);
  %    NB_real_NB_ind = sub2ind(size(grid.NBI),NBIind, INBind);
  dir_NB_ind =find(grid.NBI == -1);
  real_or_dir_NB_ind = [real_NB_ind(:); dir_NB_ind(:)];

  if model.verbose >= 29
    disp(['found ',num2str(length(real_NB_ind)),' non-boundary edges.'])
    disp(['found ',num2str(length(dir_NB_ind)),' dirichlet-boundary edges.'])
  end;
end;

% currently assuming full implicit discretization of diffusivity

%%%%%%%% 1. implicit matrix:
% has entries 
% (L_I_diff)il = 1/|T_i| sum_j (NB(i) cup NB_dir(i)) 
%                             diff*|S_ij|/d_ij            if l=i
%               -diff |S_il|/d_il                         if l is NB(i)

% evaluate all edge midpoint diffusivities at time t+dt:


if decomp_mode == 2
  diff = model.diffusivity_ptr([],model);
else
  if nargin == 4
    tmpU = repmat(U, 1, grid.nneigh);
    neiU = tmpU;
    real_nb_ind = grid.NBI > 0;
    neiU(real_nb_ind) = U(grid.NBI(real_nb_ind));

    %edgeU = 0.5 * (tmpU + neiU);
    %edgeU = sqrt(tmpU.*neiU);
    %edgeU = min(tmpU,neiU);
    %edgeU = 2*(1./tmpU + 1./neiU).^(-1);
    %diff = model.diffusivity_ptr([grid.ESX(:), grid.ESY(:)],edgeU(:), model);
    difftmp = model.diffusivity_ptr([grid.ESX(:),grid.ESY(:)],tmpU(:),model);
    diffnei = model.diffusivity_ptr([grid.ESX(:),grid.ESY(:)],neiU(:),model);

    %diff.K = 2*(difftmp.K + diffnei.K) ./ (difftmp.K + diffnei.K);
    %mymean = @harmmean;
    %mymean = @geomean;
    mymean = @mean;
    diff.K = mymean([difftmp.K, diffnei.K], 2);
    diff.epsilon = mymean([difftmp.epsilon, diffnei.epsilon]);
  else
    diff = model.diffusivity_ptr([grid.ESX(:), grid.ESY(:)],model);
  end
  % exact diffusivity on dirichlet edges
  if (~isempty(dir_NB_ind)) & (model.decomp_mode == 0)
    Xdir = grid.ECX(dir_NB_ind);
    Ydir = grid.ECY(dir_NB_ind);
    Udir = model.dirichlet_values_ptr([Xdir(:),Ydir(:)],model);
    Kdir = model.diffusivity_ptr([Xdir(:),Ydir(:)],Udir,model);
    diff.K(dir_NB_ind) = Kdir.K(:);
  end

end;


if decomp_mode == 0
  if length(diff.K)==1
    K = ones(size(grid.ESX))*diff.K;
  else
    K = reshape(diff.K,size(grid.ESX));
  end
  Ulaplace = model.laplacian_derivative_ptr([grid.CX(:),...
                                             grid.CY(:)],...
                                            U, model);
  if length(Ulaplace) == 1
    Ulaplace = ones(size(grid.CX)) * Ulaplace;
  end
  Ulaplace = repmat(Ulaplace,1,grid.nneigh);
  if ~isempty(dir_NB_ind)
    Udirlaplace = model.laplacian_derivative_ptr([Xdir,Ydir],...
                                                Udir, model);
  else
    Udirlaplace = [];
  end
  % compute products
  val = zeros(size(grid.ESX));
  val(real_or_dir_NB_ind) = K(real_or_dir_NB_ind).* ...
      grid.EL(real_or_dir_NB_ind).* ...
      [ Ulaplace(real_NB_ind); Udirlaplace ] .* ...
      (([grid.DS(real_NB_ind);1.0*grid.DS(dir_NB_ind)]).^(-1));
  Ti_inv = repmat(grid.Ainv(:),1,grid.nneigh);
  val = val .* Ti_inv;

  % diagonal values:
  L_I_diag = sparse(1:n,1:n, sum(val,2));
  % off-diagonal values:
  [i,dummy] = ind2sub(size(grid.NBI),real_NB_ind);
  L_I_offdiag = sparse(i,grid.NBI(real_NB_ind), -val(real_NB_ind),n,n);    
  L_I_diff = L_I_diag + L_I_offdiag;
  if ~isempty(NU_ind)
    L_I_diff = L_I_diff(NU_ind, :);
  end
elseif decomp_mode == 1
  % exactly identical computation as above for all components,
  % weighting as given by the diffusivity.
  L_I_diff = cell(length(diff),1);
  % temporary quantities required for all q:
  Ti_inv = repmat(grid.Ainv(:),1,grid.nneigh);
  [i,dummy] = ind2sub(size(grid.NBI),real_NB_ind);
  for q = 1:length(diff)
    K = reshape(diff{q}.K,size(grid.ESX));
    % compute products
    val = zeros(size(grid.ESX));
    val(real_or_dir_NB_ind) =  K(real_or_dir_NB_ind).* ...
        grid.EL(real_or_dir_NB_ind).* ...
        (1.0*grid.DS(real_or_dir_NB_ind)).^(-1);
    val = val .* Ti_inv;

    % diagonal values:
    L_I_diag = sparse(1:n,1:n, sum(val,2));
    % off-diagonal values:
    L_I_offdiag = sparse(i,grid.NBI(real_NB_ind), -val(real_NB_ind),n,n);    
    L_I_diff{q} = L_I_diag + L_I_offdiag;  
  end;    
else % decomp_mode==2: coefficients
     % simply forward the diffusivity sigmas
     L_I_diff = diff;
end;

%%%%%%%% 2. dirichlet-offset-vector:
% (bdir_I_diff)_i = - 1/|T_i| sum_j Ndir(i) (-d |S_ij|/d_ij u_dir(c_ij,t+dt))
if decomp_mode ==0
  if ~isempty(dir_NB_ind)
    % evaluate dirichlet values at edge midpoints at t+dt
    % Xdir, Ydir, Udir, K(dir_NB_ind) already computed above


    Kdir = K(dir_NB_ind); % diffusivity already computed above
    val2 = zeros(size(grid.ESX));
    val2(dir_NB_ind)= - grid.EL(dir_NB_ind).*...
      Udirlaplace.*((1.0*grid.DS(dir_NB_ind)).^(-1)) ...
      .* Kdir .*Udir;
    bdir_I_diff = - grid.Ainv(:).* sum(val2,2);
  else
    bdir_I_diff = zeros(grid.nelements,1);
  end
elseif decomp_mode == 1
  if ~isempty(dir_NB_ind)
    % evaluate dirichlet values at edge midpoints at t+dt
    Xdir = grid.ECX(dir_NB_ind);
    Ydir = grid.ECY(dir_NB_ind);
    Udir = model.dirichlet_values_ptr([Xdir(:),Ydir(:)],model);

    % for each combination of Udir and diffusivity component
    % perform identical computation as above
    Q_Udir = length(Udir);
    Q_diff = length(diff);
    bdir_I_diff = cell(Q_Udir * Q_diff,1);
    for q1 = 1:Q_Udir
      for q2  = 1:Q_diff
        Kdir = diff{q2}.K(dir_NB_ind); % diffusivity already computed above
        val2 = zeros(size(grid.ESX));
        val2(dir_NB_ind)= ...
          - grid.EL(dir_NB_ind).*(grid.DS(dir_NB_ind).^(-1)) ...
          .* Kdir .* Udir{q1}; 
        bdir_I_diff{(q1-1)*Q_diff+q2} = ...
          - grid.Ainv(:).* sum(val2,2);           
      end;
    end;
  else
    % still produce dummy components as the availability of dir
    % boundary cannot be detected in online phase 
    tmodel= model;
    tmodel.decomp_mode = 2;
    Udir = tmodel.dirichlet_values_ptr([],tmodel);
    % for each combination of Udir and diffusivity component
    % perform identical computation as above
    Q_Udir = length(Udir);
    Q_diff = length(diff);
    bdir_I_diff = cell(Q_Udir * Q_diff,1);
    bdir_I_diff(:) = {zeros(grid.nelements,1)};
  end;
else % decomp_mode==2 -> coefficient
  Udir = model.dirichlet_values_ptr([],model);

  % for each combination of Udir and diffusivity component
  % perform identical computation as above
  Q_Udir = length(Udir);
  Q_diff = length(diff);
  bdir_I_diff = zeros(Q_Udir * Q_diff,1);
  for q1 = 1:Q_Udir
    for q2  = 1:Q_diff
      bdir_I_diff((q1-1)*Q_diff+q2) = diff(q2) * Udir(q1);
    end;
  end;
end;

% the following cannot be detected in coefficients mode!!!
%else % no dirichlet-boundary:
%    if decomp_mode == 0
%      bdir_I_diff = zeros(n,1);
%    elseif decomp_mode == 1
%      bdir_I_diff = {};
%    else % decomp_mode == 2 -> coefficient
%      bdir_I_diff = [];
%    end;
%  end;

