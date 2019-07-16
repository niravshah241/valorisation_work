function num_flux_mat = fv_num_diff_flux(model, model_data, U, NU_ind)
%function num_flux_mat = fv_num_diff_flux(model, model_data, U, [NU_ind])
%
% function computing a numerical diffusive flux matrix for a convection
% diffusion problem. Dirichlet boundary treatment is performed, but
% neuman values are set to nan, this must be handled by calling function.
%
% In case of gradient_approximation for discretization of a
% 'gradient_tensor' the parameter NU_ind is necessary, if the gradient
% shall only be approximated on the local grid coordinates given by NU_ind.
%
% fields of num_flux:
%   epsilon: diffusion coefficient bound
%   G:  the matrix of numerical flux values across the edges.
%       boundary treatment is performed as follows:
%       in dirichlet-boundary edges, the neighbour-values are simply set
%       to this value and the flux is computed.
%
% required fields of model:
%    name_diffusive_num_flux      :  'none', 'gradient', 'gradient2'
%    'gradient':                     use edge-midpoints for diffusivity
%                                    discretization
%    'gradient2':                    use two element-cog values
%                                    for diffusivity discretization
%                                    instead of edge-midpoint value
%    'gradient_tensor':              interpolate a gradient at each edge center
%                                    and left-multiply it with a 'tensor'
%                                    matrix calculated by diffusivity_tensor.
%                                    model for diffusivity_tensor must be
%                                    given too in this case.
% plus additional fields required by dirichlet_values

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

  if isempty(NU_ind)
    NU_ind = 1:grid.nelements;
  end

  %nu_ind_length = length(NU_ind);

  num_flux_mat = [];
  num_flux_mat.epsilon = 0;
  num_flux_mat.G = zeros(size(grid.NBI));

  % determine index sets of boundary types
  real_NB_ind = find(grid.NBI>0);
%  NBIind = grid.NBI(real_NB_ind);
%  INBind = grid.INB(real_NB_ind);
%  NB_real_NB_ind = sub2ind(size(grid.NBI),NBIind,INBind);

  dir_NB_ind =find(grid.NBI == -1);
  neu_NB_ind =find(grid.NBI == -2);

  if model.verbose >= 29
    disp(['found ',num2str(length(real_NB_ind)),' non-boundary edges.'])
    disp(['found ',num2str(length(dir_NB_ind)),' dirichlet-boundary edges.'])
    disp(['found ',num2str(length(neu_NB_ind)),' neumann-boundary edges.'])
  end;

  %%%%%%%%%%%%%%%%%%%  diffusive numerical flux %%%%%%%%%%%%%%%%%%%%%%%%%%
  if isequal(model.name_diffusive_num_flux,'none')
    % set output to zero instead of nan
    num_flux_mat.G(:) = 0;
  elseif isequal(model.name_diffusive_num_flux,'gradient')
    % 1. method: evaluation of diffusivity in edge-circumcenter-intersections
    diff = model.diffusivity_ptr([grid.ESX(:),grid.ESY(:)],model);
    K = reshape(diff.K,size(grid.ESX));
    nfaces = size(grid.NBI,2);
    UU = repmat(U,1,nfaces);
    % matrix with the neighbouring values products
    UNB = nan * ones(size(UU));
    UNB(real_NB_ind) = UU(grid.NBI(real_NB_ind));

    % set real neighbour values
    num_flux_mat.G(real_NB_ind) = ...
        K(real_NB_ind).*grid.EL(real_NB_ind).*(grid.DS(real_NB_ind).^(-1).* ...
                                              (UU(real_NB_ind)- ...
                                               UNB(real_NB_ind)));

    % determine dirichlet-boundary values as required by convective
    % and diffusive flux.
    if ~isempty(dir_NB_ind)
      Xdir = grid.ESX(dir_NB_ind);
      Ydir = grid.ESY(dir_NB_ind);
      Udir = model.dirichlet_values_ptr([Xdir,Ydir],model);
      Kdir = K(dir_NB_ind);
      % determine distances circumcenter to edge-midpoint for
      % gradient approximation
      %[dir_NB_i,dir_NB_j] = ind2sub(size(UU),dir_NB_ind);
      %Ddir = sqrt((grid.ESX(dir_NB_ind)-grid.CX(dir_NB_i)).^2 + ...
        %         (grid.ESY(dir_NB_ind)-grid.CY(dir_NB_i)).^2);
      % set dirichlet neighbour values
      %num_flux_mat.G(dir_NB_ind) = ...
      %  grid.EL(dir_NB_ind).*(Ddir.^(-1).* Kdir .*(UU(dir_NB_ind)-Udir));
%      keyboard
      num_flux_mat.G(dir_NB_ind) = ...
          grid.EL(dir_NB_ind).*(grid.DS(dir_NB_ind).^(-1).* Kdir ...
                                .*(UU(dir_NB_ind)-Udir));
    end;

    % set diffusivity bound
    num_flux_mat.epsilon = diff.epsilon;

  elseif isequal(model.name_diffusive_num_flux,'gradient2')
    % 2. method: vector evaluating the diffusivity in all element
    %    circumcenters
    diff = model.diffusivity_ptr([grid.SX,grid.SY],model);

    nfaces = size(grid.NBI,2);
    UK = repmat(diff.K.*U,1,nfaces);
    % matrix with the neighbouring diffusivity products
    UKNB = nan * ones(size(UK));
    UKNB(real_NB_ind) = UK(grid.NBI(real_NB_ind));

    % set real neighbour values
    num_flux_mat.G(real_NB_ind) = ...
        grid.S(real_NB_ind).*(grid.DS(real_NB_ind).^(-1).* ...
                              (UK(real_NB_ind)- ...
                               UKNB(real_NB_ind)));

    % determine dirichlet-boundary values as required by convective
    % and diffusive flux.
    if ~isempty(dir_NB_ind)
      Xdir = grid.ESX(dir_NB_ind);
      Ydir = grid.ESY(dir_NB_ind);
      Udir = dirichlet_values(model,Xdir,Ydir);
      Kdir = diffusivity(model,Xdir, Ydir);
      % determine distances COG to edge-midpoint for gradient approximation
      %[dir_NB_i,dir_NB_j] = ind2sub(size(UK),dir_NB_ind);
      %Ddir = sqrt((grid.Sx(dir_NB_ind)-grid.CX(dir_NB_i)).^2 + ...
      %   (grid.Sy(dir_NB_ind)-grid.CY(dir_NB_i)).^2);

      % set dirichlet neighbour values
      %num_flux_mat.G(dir_NB_ind) = ...
      %  grid.S(dir_NB_ind).*(Ddir.^(-1).* (UK(dir_NB_ind)-Udir.*Kdir.K));
      num_flux_mat.G(dir_NB_ind) = ...
          grid.EL(dir_NB_ind).*(grid.DS(dir_NB_ind).^(-1).* ...
                                (UK(dir_NB_ind)-Udir.*Kdir.K));
    end;

    % set diffusivity bound
    num_flux_mat.epsilon = diff.epsilon;

  elseif isequal(model.name_diffusive_num_flux, 'gradient_tensor')
    % get values in centers of the cells adjacent to the edge an the edges
    % corner points. The corner values are computed by a 'weighted' average
    % over the four adjacent cells. Then do an lsq approximation of the
    % face the four points lie on. QUESTION: How can this be done in 3d and
    % on general grids?
    if isequal(model.gridtype, 'rectgrid')
      G_with_nans = zeros(size(num_flux_mat.G));

      for edge = 1:4
        diff   = model.diffusivity_tensor_ptr([grid.ECX(NU_ind,edge),grid.ECY(NU_ind,edge)], U, model, 1);
        diff_k = model.diffusivity_ptr([grid.ECX(NU_ind,edge),grid.ECY(NU_ind,edge)], model, U);
        % construct the diffusivity matrix D for current edge
%        diff1 = diff.K1;
%        diff2 = diff.K2;

%        [tmpgrad,bdir] = gradient_approx_matrix(model,model_data,NU_ind,edge);
        [tmpgrad] = gradient_approx(model,model_data,U,NU_ind,edge);

        % tmp2 stores the matrix vector product of 'diffusion tensor x
        % gradient'
%        tmpgrad = tmpgrad(NU_ind,:);
        vlen = size(tmpgrad,1);
        tmpgrad = reshape(tmpgrad',2*vlen,1);

%        tmpggrad = model.gravity*repmat([0;grid.EL(1,2)],vlen,1);

        tmp1    = diff_k.K * diff.K;
        tmp2    = tmp1 * tmpgrad;
%        tmp22   = tmp1 * tmpgrad2;
%        tmp2                = [ sum(diff1 .* tmpgrad, 2), ...
%                                sum(diff2 .* tmpgrad, 2) ];

%        if edge == 2 || edge == 4
%          quiver(grid.ECX(:,edge),grid.ECY(:,edge), tmpgrad(:,1), tmpgrad(:,2));
%          hold on;
%          quiver(grid.ECX(:,edge),grid.ECY(:,edge), tmp2(:,1), tmp2(:,2));
%          hold off;
%          keyboard;
%        end

        tmp2 = reshape(tmp2, 2, vlen)';

%        M = length(NU_ind);
%        offset = kron(M,NU_ind)-M;
%        indices = repmat(offset',M,1)+repmat(NU_ind,1,M);

        G_with_nans(NU_ind,edge) = sum(tmp2 .* [ grid.NX(NU_ind,edge), grid.NY(NU_ind,edge) ], 2);
%        helper = sparse([1:nu_ind_length,1:nu_ind_length], ...
%                       [2*(1:nu_ind_length)-1, 2*(1:nu_ind_length)], ...
%                        reshape( [grid.NX(NU_ind,edge), grid.NY(NU_ind,edge)], 1, 2*nu_ind_length ) );
%        G_with_nans_temp = helper * tmp2;
%        G_with_nans_temp = G_with_nans_temp * U + helper * tmp1 * bdir;
%        max(max(G_with_nans_temp - G_with_nans(NU_ind, edge)))
%        G_with_nans(NU_ind, edge) = G_with_nans_temp * U + helper * tmp1 * bdir;
      end
      G_with_nans                 = G_with_nans .* grid.EL;
      num_flux_mat.G(real_NB_ind) = G_with_nans(real_NB_ind);
      num_flux_mat.G(dir_NB_ind)  = G_with_nans(dir_NB_ind);

      num_flux_mat.epsilon        = max(max(num_flux_mat.G));
    else
      error(['gradient_tensor is not implemented for specified gridtype ', ...
              model.gridtype]);
    end
  else
    error ('diffusive numerical flux unknown');
  end;

% vim: set sw=2 et:
%| \docupdate 
