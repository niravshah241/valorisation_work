function [L_I_diff, bdir_I_diff] = fv_operators_diff_implicit_gradient_tensor(...
    model,model_data,NU_ind)
%function [L_I_diff, bdir_I_diff] = ...
%    fv_operators_diff_implicit_gradient_tensor(model,model_data[,NU_ind])
%
% function computing the implicit diffusion contribution to the time evolution
% matrices for a finite volume time step  L_I * Unew = L_E * U + b
% this function is called by fv_operators_implicit_explicit.
%
% the *_implicit functions perform a dt increase in model
% themselves before evaluating data functions
%
% Result is a sparse matrix L_I_diff and an offset vector
% bdir_I_diff, the latter containing dirichlet value contributions
%
% as numerical flux simple gradient approximation is used
% 
% Function supports affine decomposition, i.e. different operation modes
% guided by optional field decomp_mode in model. See also the 
% contents.txt for general explanation


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


old_t = model.t;
model.t=model.t+model.dt;
  
if decomp_mode == 2 % coefficients
  % simply forward the diffusivity sigmas
  L_I_diff    = [];
  bdir_I_diff = [];
  for edge = 1:4
    [tmpgrad,bdir] = gradient_approx_matrix(model, [], [], edge);
    diff_k         = model.diffusivity_ptr([], model);
    L_I_temp       = kron(tmpgrad, diff_k);
    b_I_temp       = kron(bdir, diff_k);
    L_I_diff       = [ L_I_diff, L_I_temp ];
    bdir_I_diff    = [ bdir_I_diff, b_I_temp ];
  end
else
  grid = [];
  if ~isempty(model_data)
    grid = model_data.grid;
  end

  %nu_ind_set = false;

  if nargin < 3 || isempty(NU_ind)
    NU_ind = 1:grid.nelements;
%  else
    %nu_ind_set = true;
  end

  nu_ind_length = length(NU_ind);
  n = grid.nelements;
  % determine index sets of boundary types
  real_NB_ind = find(grid.NBI>0);
  % NBIind = grid.NBI(real_NB_ind);
  % INBind = grid.INB(real_NB_ind);
  %    NB_real_NB_ind = sub2ind(size(grid.NBI),NBIind, INBind);
  dir_NB_ind =find(grid.NBI == -1);
  % real_or_dir_NB_ind = [real_NB_ind(:); dir_NB_ind(:)];

  if model.verbose >= 29
    disp(['found ',num2str(length(real_NB_ind)),' non-boundary edges.'])
    disp(['found ',num2str(length(dir_NB_ind)),' dirichlet-boundary edges.'])
  end;

  if ~isequal(model.gridtype, 'rectgrid')
    error(['gradient_tensor is not implemented for specified gridtype ', ...
      model.gridtype]);
  end

  U = [];

  if decomp_mode == 0 % complete mode
    tmp_flux_mat = sparse(nu_ind_length, n);
    bdir_I_diff  = zeros(nu_ind_length, 1);
    for edge = 1:4
      [tmpgrad,bdir] = gradient_approx_matrix(model, model_data, NU_ind, edge);
      diff_k         = model.diffusivity_ptr([grid.ECX(NU_ind,edge),grid.ECY(NU_ind,edge)], model);
      diff           = ...
        model.diffusivity_tensor_ptr([grid.ECX(NU_ind,edge), grid.ECY(NU_ind,edge)], ...
                                     U, model, 3);
 
%      if nu_ind_set && model.debug
%        disp(['any wrong: ', num2str(any(bdir2(reshape([2*NU_ind-1,2*NU_ind]',2*length(NU_ind),1))-bdir))]);
%  
%        if any(bdir2(reshape([2*NU_ind-1,2*NU_ind]',2*length(NU_ind),1))-bdir)
%          [bdir, bdir2(reshape([2*NU_ind-1,2*NU_ind]',2*length(NU_ind),1))]
%          keyboard;
%        end
%      end

      % this works only for scalar diffusivity -> use only first value of
      % diffusivity vector (assumed to be homogeneous)
      if model.debug
        if any(repmat(diff_k.K(1),length(diff_k.K),1) ~= diff_k.K)
          error(['Error: model.diffusivity_ptr does not return a homogeneous vector!\n',...
                 'This is invalid inside fv_operators_diff_implicit_gradient_tensor.']);
        end
      end
      tmp1   = diff_k.K(1) * diff.K;

      helper = sparse([1:nu_ind_length,1:nu_ind_length], ...
                      [2*(1:nu_ind_length)-1, 2*(1:nu_ind_length)], ...
                      reshape( [grid.NX(NU_ind,edge), grid.NY(NU_ind,edge)], ...
                      1, 2*nu_ind_length ) );
  
      tmp2 = helper * tmp1;
  
      ELL = spdiags(grid.EL(NU_ind,edge).*grid.Ainv(NU_ind,1),0,nu_ind_length,nu_ind_length);
      tmp_flux_mat = tmp_flux_mat + ELL * tmp2 * tmpgrad;
      bdir_I_diff  = bdir_I_diff - tmp2 * bdir .* grid.EL(NU_ind,edge).*grid.Ainv(NU_ind,1);
    end
    L_I_diff = - tmp_flux_mat;
  elseif decomp_mode == 1 % (components mode)
    L_I_diff    = cell(0,1);
    bdir_I_diff = cell(0,1);
    for edge = 1:4
      model.decomp_mode = 0;
      diff           = model.diffusivity_tensor_ptr([grid.ECX(NU_ind,edge), ...
                                                     grid.ECY(NU_ind,edge)],[],model, 3);
  
      model.decomp_mode = 1;
      [tmpgrad,bdir] = gradient_approx_matrix(model, model_data, NU_ind, edge);
      diff_k         = model.diffusivity_ptr([grid.ECX(NU_ind,edge),grid.ECY(NU_ind,edge)], model);
  
      L_I_temp = cell(length(diff_k) * length(tmpgrad), 1);
      b_I_temp = cell(length(diff_k) * length(bdir), 1);
  
      helper = sparse([1:nu_ind_length,1:nu_ind_length], ...
                      [2*(1:nu_ind_length)-1, 2*(1:nu_ind_length)], ...
                      reshape( [grid.NX(NU_ind,edge), grid.NY(NU_ind,edge)], ...
                      1, 2*nu_ind_length ) );
  
      for q=1:length(diff_k)
        tmp1 = spdiags(repmat(diff_k{q}.K,2,1),0,size(diff.K,1),size(diff.K,2)) * diff.K;
        tmp2 = helper * tmp1;
        for r=1:length(tmpgrad)
          %          t = diag(grid.EL(:,edge)) * tmp2 * tmpgrad{r};
          ELL = spdiags(grid.EL(NU_ind,edge).*grid.Ainv(NU_ind,1),0,nu_ind_length,nu_ind_length);
          L_I_temp{ (q-1)*length(tmpgrad) + r } = ELL * tmp2 * tmpgrad{r};
        end
        for r=1:length(bdir)
          b_I_temp{ (q-1)*length(bdir) + r } = tmp2 * bdir{r} .* grid.EL(NU_ind,edge).*grid.Ainv(NU_ind,1);
        end
      end
      L_I_diff    = [ L_I_diff; L_I_temp ];
      bdir_I_diff = [ bdir_I_diff; b_I_temp ];
    end
  else
    error(['decomp_mode number ', model.decomp_mode, ' is unknown.']);
  end

end

model.t = old_t;

