function num_flux = fv_num_diff_flux_gradient_tensor(model,model_data,U, NU_ind)
%function num_flux = fv_num_diff_flux_gradient_tensor(model,model_data,U, NU_ind)
% computes a numerical diffusive flux for a diffusion problem including a tensor
%
% This function computes a numerical diffusive flux matrix for a diffusion
% problem. Dirichlet boundary treatment is performed, but neuman values are set
% to 'NaN'. These must be handled by the calling function.
%
% For discretization of a 'gradient_tensor' the parameter NU_ind is
% necessary, if the gradient shall only be approximated on the local
% grid coordinates given by NU_ind.
%
% The analytical term inspiring this flux looks like
% ` - \nabla \cdot \left( d(x) K(x,u) \nabla u \right),`
% with a function `d: \mathbb{R} \to \mathbb{R}` and a tensor `K: R \times R \to \text{Lin}(R^2, R^2)`.
%
% The implemented flux then can be written for inner boundaries as
%    ``g_{ij}(u,v) = d(x_{ij}) \tilde{K}(x,u,v)
%    \tilde{\nabla}(u,v) |e_{ij}|``
% and for Dirichlet boundaries as
%    ``g_{dir,ij}(u) = d(x_{ij}) \tilde{K}_{dir}(x,u) \tilde{\nabla}(u) |e_{ij}|``
% where the discretizations ``\tilde{K} and \tilde{\nabla}``
% of K
% respectively `\nabla` are given by the function pointer
% 'model.diffusivity_tensor_ptr' respectively the function gradient_approx().
% Note, that the latter is only available for rectangular grids.
%
% required fields of model:
%   diffusivity_ptr        : function pointer to diffusivity function `d`
%   diffusivity_tensor_ptr : function pointer to a tensor implementation for `K`
%
% generated fields of num_flux:
%   epsilon: diffusion coefficient bound
%   G:  the matrix of numerical flux values across the edges.
%       boundary treatment is performed as follows:
%       in dirichlet-boundary edges, the neighbour-values are simply set
%       to this value and the flux is computed.
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

if isempty(NU_ind)
  NU_ind = 1:grid.nelements;
end

%nu_ind_length = length(NU_ind);

num_flux = [];
num_flux.epsilon = 0;
num_flux.G = zeros(size(grid.NBI));

% determine index sets of boundary types
real_NB_ind = find(grid.NBI>0);
dir_NB_ind =find(grid.NBI == -1);
neu_NB_ind =find(grid.NBI == -2);

if model.verbose >= 29
  disp(['found ',num2str(length(real_NB_ind)),' non-boundary edges.'])
  disp(['found ',num2str(length(dir_NB_ind)),' dirichlet-boundary edges.'])
  disp(['found ',num2str(length(neu_NB_ind)),' neumann-boundary edges.'])
end

% get values in centers of the cells adjacent to the edge an the edges
% corner points. The corner values are computed by a 'weighted' average
% over the four adjacent cells. Then do an lsq approximation of the
% face the four points lie on. QUESTION: How can this be done in 3d and
% on general grids?
if model.debug && ~isequal(model.gridtype, 'rectgrid')
  error('gradient_tensor only implemented for rectangular grids');
end

G_with_nans = zeros(size(num_flux.G));

for edge = 1:4
  diff = model.diffusivity_tensor_ptr([grid.ECX(NU_ind,edge), ...
                                       grid.ECY(NU_ind,edge)], ...
                                      U, model, 1);
  model.U = U(NU_ind);
  diff_k = model.diffusivity_ptr([grid.ECX(NU_ind,edge),...
                                  grid.ECY(NU_ind,edge)], model);

  % construct the diffusivity matrix D for current edge
  %[tmpgrad,bdir] = gradient_approx_matrix(model,model_data,...
  %                                        NU_ind,edge);
  [tmpgrad] = gradient_approx(model,model_data,U,NU_ind,edge);

  % tmp2 stores the matrix vector product of 'diffusion tensor x
  % gradient'
  vlen = size(tmpgrad,1);
  tmpgrad = reshape(tmpgrad',2*vlen,1);

  % compute on each intersection e = e_{i,edge} the diffusion tensor 
  %   D_e = d_e * [ T_{e11}, T_{e12};
  %                 T_{e21}, T_{e22} ] ...
  tmp1    = spdiags(repmat(diff_k.K,2,1),0,size(diff.K,1),size(diff.K,2)) * diff.K;
  % ... and multiply it with approximated gradient over the intersection 
  %    g_e \circ \nabla U | _e:
  %    h_e = D_e * g_e
  % The result is
  %    tmp2 = [ h_{e_{1,edge}}; h_{e_{2,edge}}; ... ; h_{e_{H,edge}} ].
  tmp2    = tmp1 * tmpgrad;

%   if edge == 2 || edge == 4
%     quiver(grid.ECX(:,edge),grid.ECY(:,edge), tmpgrad(:,1), tmpgrad(:,2));
%     hold on;
%     quiver(grid.ECX(:,edge),grid.ECY(:,edge), tmp2(:,1), tmp2(:,2));
%     hold off;
%     keyboard;
%   end

  tmp2 = reshape(tmp2, 2, vlen)';

  % M = length(NU_ind);
  % offset = kron(M,NU_ind)-M;
  % indices = repmat(offset',M,1)+repmat(NU_ind,1,M);

  G_with_nans(NU_ind,edge) = sum(tmp2 .* [ grid.NX(NU_ind,edge),...
                                           grid.NY(NU_ind,edge) ], 2);
end
G_with_nans             = - G_with_nans .* grid.EL;
num_flux.G(real_NB_ind) = G_with_nans(real_NB_ind);
num_flux.G(dir_NB_ind)  = G_with_nans(dir_NB_ind);

num_flux.epsilon        = max(max(num_flux.G));

