function num_flux_mat = fv_num_diff_flux_gradient(model, model_data, U, NU_ind)
%function num_flux_mat = fv_num_diff_flux_gradient(model, model_data, U, [NU_ind])
% computes a numerical diffusive flux for a diffusion problem
%
% This function computes a numerical diffusive flux matrix for a diffusion
% problem. Dirichlet boundary treatment is performed, but neuman values are set
% to 'NaN'. These must be handled by the calling function.
%
% The analytical term inspiring this flux looks like
% `` \nabla \cdot \left(d(u) \nabla l(u) \right). ``
% Here, the functions `d(u)` and `l(u)` differ in the numerical implementation.
% The latter is evaluated in cell centers during flux computation, whereas the
% diffusion functional `d(u)` is evaluated on the cell edges by averaging the
% solutions value there.
%
% The implemented flux then can be written for inner boundaries as
%   `` g_{ij}(u,v) = d(\{u,v\}) \frac{l(v) - l(u)}{|d_{ij}|} |e_{ij}| ``
% and for Dirichlet boundaries as
%   `` g_{dir,ij}(u) = d(u_{dir}(x_{ij}))
%          \frac{l(u_{dir}(x_{ij})) - l(u)}{|d_{ij}|} |e_{ij}|, ``
% where `\{u,v\}` refers to the average of `u` and `v`.
%
% required fields of model:
%   diffusivity_ptr:  function pointer to diffusivity function `d`
%   laplacian_ptr:    function pointer to nonlinear function `l`
%
% generated fields of num_flux_mat:
%   epsilon: diffusion coefficient bound
%   G:  the matrix of numerical flux values across the edges.
%       boundary treatment is performed as follows:
%       in dirichlet-boundary edges, the neighbour-values are simply set
%       to this value and the flux is computed.
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

if nargin == 3 || isempty(NU_ind)
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

% 1. method: evaluation of diffusivity in edge-circumcenter-intersections

tmpU = repmat(U, 1, grid.nneigh);
neiU = tmpU;
real_nb_ind = grid.NBI > 0;
neiU(real_nb_ind) = U(grid.NBI(real_nb_ind));
%edgeU = 0.5 * (tmpU + neiU);
%edgeU = sqrt(tmpU.*neiU);
%edgeU = min(tmpU, neiU);
%edgeU = 2*(tmpU.*neiU)./(tmpU+neiU);

model.U = tmpU(:);
difftmp = model.diffusivity_ptr([grid.ESX(:),grid.ESY(:)],model);
model.U = neiU(:);
diffnei = model.diffusivity_ptr([grid.ESX(:),grid.ESY(:)],model);
%mymean = @harmmean;
%mymean = @geomean;
mymean = @mean;
diff.K = mymean([difftmp.K, diffnei.K],2);
diff.epsilon = mymean([difftmp.epsilon, diffnei.epsilon]);
if length(diff.K) == 1
  K = ones(size(grid.ESX))*diff.K;
else
  K = reshape(diff.K,size(grid.ESX));
end
nfaces = size(grid.NBI,2);

UU = model.laplacian_ptr([grid.CX(:), grid.CY(:)], U(:), model);

UU = repmat(UU,1,nfaces);
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
  Xdir  = grid.ESX(dir_NB_ind);
  Ydir  = grid.ESY(dir_NB_ind);
  Udir  = model.dirichlet_values_ptr([Xdir,Ydir],model);
  Uldir = model.laplacian_ptr([Xdir,Ydir],Udir,model);
  model.U = Udir;
  Kdir  = model.diffusivity_ptr([Xdir,Ydir], model);
  Kdir  = Kdir.K;
  % determine distances circumcenter to edge-midpoint for
  % gradient approximation
  %[dir_NB_i,dir_NB_j] = ind2sub(size(UU),dir_NB_ind);
  %Ddir = sqrt((grid.ESX(dir_NB_ind)-grid.CX(dir_NB_i)).^2 + ...
    %         (grid.ESY(dir_NB_ind)-grid.CY(dir_NB_i)).^2);
  % set dirichlet neighbour values
  %num_flux_mat.G(dir_NB_ind) = ...
  %  grid.EL(dir_NB_ind).*(Ddir.^(-1).* Kdir .*(UU(dir_NB_ind)-Udir));
  num_flux_mat.G(dir_NB_ind) =  ...
      grid.EL(dir_NB_ind).*((1.0*grid.DS(dir_NB_ind)).^(-1).* Kdir ...
                            .*(UU(dir_NB_ind)-Uldir));
end

% set diffusivity bound
num_flux_mat.epsilon = diff.epsilon;

end

% vim: set sw=2 et:
%| \docupdate 
