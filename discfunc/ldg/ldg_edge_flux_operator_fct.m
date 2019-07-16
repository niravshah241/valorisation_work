function L = ldg_edge_flux_operator_fct(llcoord,grid,params)
%function L = ldg_edge_flux_operator_fct(llcoord,grid,params)
%
% function returning the integrand of the edge flux integral
% evaluated in 1D local coordinates llcoord along the reference
% interval
%
% in 'components' mode, a 5 dimensional dense matrix is generated
% L(i,j,e,f,q) = v^q (lcoord(e,f,llcoord)) * n_(e,f) *
%                   phi_e_i(lcoord(e,f,llcoord))
%                   phi_e_j (lcoord(e,f,llcoord))
%                 * |f|
% i.e. i and j local base function indices, e the element index, f
% the face/edge index and q the component number. n_(e,f) the outer
% normal of face f wrt element e, |f| the length of edge f and 
% lcoord(e,f,llcoord) the 2d local coordinates of the point in the
% reference triangle. 
%
% in 'complete' mode, the above Q = Q_v components are weighted
% with Theta_v^q(mu), i.e. identical coefficients as velocity.
% result is a 4 dimensional matrix
% L(i,j,e,f) = v (lcoord) * n_(e,f)* phi_e_i(lcoord) phi_e_j (lcoord)
%                 * |f|
%
% intervalquadrature of this function over the edges should give the
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

disp(['to be adjusted!!! ']);
keyboard;

decomp_mode = params.decomp_mode

if params.dimrange~=1
  error('currently only scalar valued implementation supported');
end;
if (decomp_mode == 2)
  L = params.velocity_coefficients_fct(params);
else
  % loop over edge types
  %determine 2d lcoord for the 3 edge types
  lcoords = llocal2local(grid,1:3,llcoord);

  if decomp_mode == 1
    L = zeros(nbasefct, nbasefct, nelements, 3, Q_v);
  else % complete
    L = zeros(nbasefct, nbasefct, nelements, 3);
  end;
  
  for f = 1:3
    % evaluate reference basis in 3 lcoords
    Phi_hats = ldg_evaluate_basis(lcoord(:,f),params);
    nbasefct = size(Phi_hats,2);
    nelements = grid.nelements;
    % evaluate velocity in lcoord of all elements
    v = params.velocity_ptr(...
	local2global(grid,1:nelements,lcoords(:,f),params), params);
    
    if (decomp_mode == 1)

      % multiply with outer unit normal
      Q_v = length(v);
      for q=1:Q_v
	vn = v{q}(1,:).*grid.NX(:,f) + v{q}(2,:).*grid.NY(:,f);
	% expand vn{q} over i range and j range
        vv = repmat(reshape(vn,[1,1,nelements]),[nbasefct,nbasefct,1]);
	% expand phi_i over j range and vn range
        pphi_i = repmat(reshape(Phi_hats),[nbasefunc,1,1],...
			[1,nbasefunc,nelements]);
	% expand phi_j over i range and vn range
        pphi_j = repmat(reshape(Phi_hats),[1,nbasefunc,1],...
			[nbasefunc,1,nelements]);
	% multiply and distribute into output matrix
        L(:,:,:,f,q) = reshape(vv.*pphi_i.*pphi_j,...
			       [numbasefct,numbasefct,nelements]);
      end;    
      
    else % decomp_mode = 0 === "complete"
      
      % multiply with outer unit normal
      vn = v(1,:).*grid.NX(:,f) + v(2,:).*grid.NY(:,f);
      % expand vn{q} over i range and j range
      vv = repmat(reshape(vn,[1,1,nelements]),[nbasefct,nbasefct,1]);
      % expand phi_i over j range and vn range
      pphi_i = repmat(reshape(Phi_hats),[nbasefunc,1,1],...
		      [1,nbasefunc,nelements]);
      % expand phi_j over i range and vn range
      pphi_j = repmat(reshape(Phi_hats),[1,nbasefunc,1],...
		      [nbasefunc,1,nelements]);
      % multiply and distribute into output matrix
      L(:,:,:,f) = reshape(vv.*pphi_i.*pphi_j,...
			   [numbasefct,numbasefct,nelements]);
      
    end;  % else
  end; % for f=1:3 loop
  
end;%| \docupdate 
