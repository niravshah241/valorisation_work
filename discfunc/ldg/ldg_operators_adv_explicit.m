[L_E, b_dir]  = ldg_operators_adv_explicit(model,model_data)
%[L_E, b_dir]  = ldg_operators_adv_explicit(model,model_data)
%
% function computing the explicit operator matrices for an
% advection problem with ldg discretization. As numerical flux,
% simple upwinding is used currently. Can be more refined by other
% numerical fluxes.

% Bernard Haasdonk 30.8.2009

disp('to be adjusted!!');
keyboard;

decomp_mode = model.decomp_mode;

if decomp_mode ==2
  % no not forget the minus - theta^q(v) for the analytical 
  % flux components!
  v_coeff = params.velocity_coefficients_fct(params);  
  num_flux_matrix_coeff = ldg_edge_num_flux_matrix([],params);
  L_E = [v_coeff(:);num_flux_matrix_coeff(:)];
  b_dir = ...
  error('to_be_implemented')
  
elseif decomp_mode == 1
  % compute inverted mass matrix, which is multiplied
  % to all L_E and b_dir components later
  params.dimrange = model.dimrange;
  params.pdeg = model.pdeg;  
  Minv = ldg_inv_mass_matrix(params,model_data.grid,model.qdeg);
  nbasefct =  size(Minv,1);
  nelements = model_data.grid.nelements;
  nfaces = 3; % triagrid
  
  % initialize outputs
  model.decomp_mode = 2; % coefficients
  Q_v = length(model.velocity_ptr([],model));
  Q_num_flux = length(ldg_adv_num_flux_matrix(model));
  Q_dir = length(model.dirichlet_values_ptr([],model));
  model.decomp_mode = 1; % components again

  L_E = cell(1,Q_v+Q_num_flux);
  b_dir = cell{1,Q_num_flux*Q_dir};
  
  % 1. components of L_E: element flux integral
  % all components of v
  % entries of L_E_q: 
  % (L_E_q)_(e,i)(e,j) = int_e v_q grad phi_e,i phi_e,j 
  
  % Integrate 4D-matrix function, but only once evaluating basis
  L = triaquadrature(model.element_flux_qdeg,...
		     @ldg_element_flux_operator_fct,...
		     model_data.grid,params);
  % multiply all integrals by Det(DF) in order to have element
  % integral instead of reference element integral
  DetDF = repmat(reshape(grid.A*2,[1,1,nelements,1]),...
		 nbasefct,nbasefct,1,Q_v);
  L = L.*DetDF;
  
  % distribute matrix layers into L_E sparse matrix components
  for q = 1:Q_v
    L_E{q} = spblkdiag(L(:,:,:,q));
  end;
  
  % 2. components of L_E: numerical flux
  % all components of ldg_adv_num_flux_matrix 
  % (== comp of v for upwind flux)

  edge_flux_mat = intervalquadrature(model.edge_flux_qdeg,...
				     @ldg_edge_flux_operator_fct,...
				     model_data.grid,params);
  %multiply with edgelengths in order to obtain real edge integral  
  EL = repmat(reshape(model_data.grid.EL,[1,1,nelements,nfaces,1]),...
	      [nbasefct,nbasefct,1,1,Q_v]);
  edge_flux_mat = edge_flux_mat .* EL;

  % evaluate numerical flux matrix components from edge-flux-matrix
  edge_num_flux_mat = ldg_edge_num_flux_matrix(edge_flux_mat,model);

  ...
  
  % distribute into further L_E components  
  for q = 1:Q_num_flux
    L_E{q+Q_v} = sparse(size(Minv));    
    ... 
    
  end;
  
  % 3. components of b_dir
  % all combinations of bdir and num_flux_matrix components
  for q1 = 1:Q_num_flux
    for q2 = 1:Q_dir
      b_dir{(q1-1)*Q_dir+q2} = zeros(size(Minv,1),1);
    end;
  end;
  ... % efficient integration over all elements at once but only
      % once evaluating basis functions
  
  % invert, such that explicit time discretization can be used
  for q = 1:length(L_E)
    L_E{q} = Minv * L_E{q};
  end;
  for q = 1:length(b_dir)
    b_dir{q} = Minv * b_dir{q};
  end;

else % "complete" == decomp_mode = 0 
  % perform above for complete velocity field
  % should not be a simple linear combination, as a new computation
  % will be very surely faster. 
  disp('complete evaluation should be accelerated in ldg_operators_adv_explicit');
  model.decomp_mode = 2;
  [L_E_coeff, b_dir_coeff] = ldg_operators_adv_explicit(model,model_data);
  model.decomp_mode = 1;
  [L_E_comp, b_dir_comp] = ldg_operators_adv_explicit(model,model_data);
  model.decomp_mode = 0;
  L_E = lincomb_sequence(L_E_comp,L_E_coeff);
  b_dir = lincomb_sequence(b_dir_comp,b_dir_coeff);  
end;

%| \docupdate 
