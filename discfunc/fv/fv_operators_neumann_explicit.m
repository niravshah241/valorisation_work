function [L_E_neu, b_E_neu] = fv_operators_neumann_explicit(model,model_data,U,NU_ind)
%function [L_E_neu, b_E_neu] = fv_operators_neumann_explicit(model,model_data,[U,NU_ind])
% computes a neumann contribution matrix for finite volume time evolution operators,
% <b> or their Frechet derivative </b>
%
% This function computes the neumann boundary contribution operator
% `L_{\text{neu}}` and a corresponding vector ` b_{\text{neu}}` that can be used
% by fv_operators_implicit_explicit() to build evolution matrices for a finite
% volume time step `L_I U^{k+1} = L_E U^k + b_E + b_I`. This operator
% contribution must be activated via the flag
% 'model.operators_neumann_implicit'
%
% The analytical terms inspiring this operator look like
% ` f(u) \cdot n = b_{\text{neu}}(u) \quad
%   \text{on }\partial \Omega_{\text{neu}},`
% where `f` can be some flux function, e.g. chosen via the
% 'model.conv_flux_ptr' and `b_{\text{neu}}(u)` is the neumann boundary
% condition selected via 'model.neumann_values_ptr'.
%
% See also: fv_operators_conv_explicit_engquist_osher(),
% fv_operators_diff_implicit_gradient(),
% fv_operators_diff_implicit_gradient_tensor()
%
% Required fields of model:
%  neumann_values_ptr :
%
% Return values:
%  L_E_neu : sparse matrix `L_{\text{neu}}`
%  b_E_neu : offset vector `b_{\text{neu}}`

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

if (decomp_mode==2)
  % only get and forward sigmas:
  bneu0 = model.neumann_values_ptr([],[],[],model);
  L_E_neu = bneu0;
  b_E_neu = bneu0;

else % decomp_mode == 0 or 1

  n = grid.nelements;
  % determine index sets of boundary types
  real_NB_ind = find(grid.NBI>0);
  %  dir_NB_ind =find(grid.NBI == -1);
  neu_NB_ind =find(grid.NBI == -2);

  if model.verbose >= 10
    disp(['found ',num2str(length(real_NB_ind)),' non-boundary edges.'])
    disp(['found ',num2str(length(neu_NB_ind)),' neumann-boundary edges.'])
  end;

  if ~isempty(neu_NB_ind)
    % determine linearization of neumann boundary values
    % FNneu = aneu * u + bneu
    Xneu = grid.ECX(neu_NB_ind);
    Yneu = grid.ECY(neu_NB_ind);

    if model.flux_linear
      Uneu = zeros(length(Xneu),1);
      bneu0 = model.neumann_values_ptr([Xneu(:),Yneu(:)],...
                                       Uneu,...
                                       [grid.NX(neu_NB_ind),grid.NY(neu_NB_ind)],...
                                       model);

      Uneu = ones(length(Xneu),1);
      bneu1 = model.neumann_values_ptr([Xneu(:),Yneu(:)],...
                                       Uneu,...
                                       [grid.NX(neu_NB_ind),grid.NY(neu_NB_ind)],...
                                       model);
      if decomp_mode == 0
        %%%%%%%% explicit matrix:
        % has entries
        % (L_E_neu)il = 1/|Ti| *
        %          (sum j NBneu(i) |S_ij| a_neu_ij       for l=i
        %                                                       0       else

        aneu = bneu1-bneu0;

        val = zeros(size(grid.ESX));
        val(neu_NB_ind)= grid.EL(neu_NB_ind).* aneu;
        L_E_neu = sparse(1:n,1:n,  grid.Ainv(:).* sum(val,2));

        % check nonnegativity
        if ~isempty(find(diag(L_E_neu)<0,1))
          disp('warning: neumann-boundary has inflow velocity !!!! ');
        end;
        %    keyboard;

        %%%%%%%% offset-vector:
        % (b_E_neu)_i = - 1/|T_i| * ...
        %    sum_j NBneu(i) |S_ij| b_E_neu_ij
        val(neu_NB_ind) = grid.EL(neu_NB_ind).* bneu0;
        b_E_neu         = - grid.Ainv(:) .* sum(val,2);
      else % decomp_mode == 1 % identical computation for all components
        Q_bneu = length(bneu0);
        L_E_neu = cell(Q_bneu,1);
        b_E_neu = cell(Q_bneu,1);
        for q=1:Q_bneu
          aneu = bneu1{q}-bneu0{q};
          val = zeros(size(grid.ESX));
          val(neu_NB_ind)= grid.EL(neu_NB_ind).* aneu;
          L_E_neu{q} = sparse(1:n,1:n,  grid.Ainv(:).* sum(val,2));
          val(neu_NB_ind)= grid.EL(neu_NB_ind).* bneu0{q};
          b_E_neu{q}= - grid.Ainv(:) .* sum(val,2);
        end;
      end; % decomp more select
    else  % compute frechet derivative ( non-linear flux )
      old_conv_flux_ptr = model.conv_flux_ptr;
      model.conv_flux_ptr = model.conv_flux_derivative_ptr;

      bneu = model.neumann_values_ptr([Xneu(:),Yneu(:)], U,...
                                      [grid.NX(neu_NB_ind), grid.NY(neu_NB_ind)],...
                                      model);
      model.conv_flux_ptr = old_conv_flux_ptr;
      if decomp_mode ~= 0
        error('only decomp mode == 0 is allowed!')
      end
      val = zeros(size(grid.ESX));
      val(neu_NB_ind)= grid.EL(neu_NB_ind).* bneu;
      L_E_neu = sparse(1:n,1:n,  grid.Ainv(:).* sum(val,2));
      b_E_neu = 0;

      % check nonnegativity
      if ~isempty(find(diag(L_E_neu)<0,1))
        disp('warning: neumann-boundary has inflow velocity !!!! ');
      end;
    end
  else % neumann boundary empty
    if decomp_mode == 0
      L_E_neu = sparse(n,n);
      b_E_neu = zeros(n,1);
    else % decomp_mode == 1
         %      L_E_neu = {}; %{sparse(n,n)};
         %      b_E_neu = {}; % {zeros(n,1)};
         % cannot determine in mode==2, whether neumann boundary is empty. So
         % at least one component must be provided.
         L_E_neu = {sparse(n,n)};
         b_E_neu = {zeros(n,1)};
    end; % neumann boundary nonempty
  end; % neumann-boundary empty select

end; % decomp_mode ==0 or 1


