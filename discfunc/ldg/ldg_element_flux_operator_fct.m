function L = ldg_element_flux_operator_fct(lcoord,grid,params)
%function L = ldg_element_flux_operator_fct(lcoord,grid,params)
%
% function returning the integrand of the element flux integral
% evaluated in local coordinates lcoord.
%
% in 'components' mode, a 4 dimensional dense matrix is generated
% L(i,j,e,q) = v^q (lcoord) * JIT* grad phi_e_i(lcoord) phi_e_j (lcoord)
%                 * |Det(DF)|
% i.e. i and j local base function indices, e the element index and 
% q the component number.
%
% in 'complete' mode, the above Q = Q_v components are weighted
% with Theta_v^q(mu), i.e. identical coefficients as velocity.
% result is a 3 dimensional matrix
% L(i,j,e) = v (lcoord) * JIT* grad phi_e_i(lcoord) phi_e_j (lcoord)
%                 * |Det(DF)|
%
% triaquad of this function over the element should give the
% desired result used in ldg_operaotors_adv_explicit

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


% Bernard Haasdonk 31.8.2009

disp([' to be adjusted!!  ']);
keyboard;

decomp_mode = params.decomp_mode

if params.dimrange~=1
  error('currently only scalar valued implementation supported');
end;
if (decomp_mode == 2)
  L = params.velocity_coefficients_fct(params);
else
  % evaluate reference basis in lcoord
  Phi_hats = ldg_evaluate_basis(lcoord,params);
  nbasefct = size(Phi_hats,2);
  nelements = grid.nelements;
  Phi_hats_derivative = ldg_evaluate_scalar_basis_derivative(lcoord,params);
  % extrude Phi_hats_derivative to all elements
  %Phi_hats_derivative = reshape([Phi_hats_derivative{:}],2, ...
  %				nbasefct);
  Phi_hats_derivative = repmat(Phi_hats_derivative,1,grid.nelements);
  JIT = repmat(grid.JIT(:)',nbasefct,1);
  JIT = reshape(JIT,nbasefct*nelements,2,2);
    
  % multiply Phi_hats_derivative with JIT (JacobianInverseTransposed)
  Phi_derivative_x = JIT(:,1,1)*Phi_hats_derivative(1,:) ...
      + JIT(:,1,2)*Phi_hats_derivative(2,:);
  Phi_derivative_y = JIT(:,2,1)*Phi_hats_derivative(1,:) ...
      + JIT(:,2,2)*Phi_hats_derivative(2,:);
  Phi_derivative = [Phi_derivative_x; Phi_derivative_y];
  
  % evaluate velocity in lcoord
  v = params.velocity_ptr(...
      local2global(grid,1:grid.nelements,lcoord,params), params);

  if (decomp_mode == 1)
    Q_v = length(v);
    L = zeros(nbasefct, nbasefct, nelements, Q_v);
    for q=1:Q_v
      % extrude v_q, such that simple componentwise multiplication
      % with others is possible
      vq = repmat(v{q}(:),1,nbasefct);
      vq = reshape(vq,2,nbasefct*nelements);
      % sort results into L: First v * grad Phi_e_i in identical columns
      L0 = v(1,:).*Phi_derivative(1,:)+ v(2,:).*Phi_derivative(2,:);
      L0 = reshape(L0,[nbasefct,1,nelements]);
      Lq = repmat(L0,1,nbasefct,1);
      % then multiply with phi_e_j   
      for j = 1: nbasefct
	L(:,j,:,q) = Lq(:,j,:).*Phi_hats(:,j);
      end;    
    end;    
    
  else % decomp_mode = 0 === "complete"
  
    % extrude v_q, such that simple componentwise multiplication
    % with others is possible
    v = repmat(v(:),1,nbasefct);
    v = reshape(v,2,nbasefct*nelements);
    % sort results into L: First v * grad Phi_e_i in identical columns
    L0 = v(1,:).*Phi_derivative(1,:)+ v(2,:).*Phi_derivative(2,:);
    L0 = reshape(L0,[nbasefct,1,nelements]);
    L = repmat(L0,1,nbasefct,1);
    % then multiply with phi_e_j   
    for j = 1: nbasefct
      L(:,j,:) = L(:,j,:).*Phi_hats(:,j);
    end;    
  
  end;  
end;%| \docupdate 
