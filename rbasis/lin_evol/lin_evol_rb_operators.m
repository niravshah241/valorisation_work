function [LL_I, LL_E, bb, K_II, K_IE, K_EE, m_I, m_E, m] = ...
    lin_evol_rb_operators(model, detailed_data)
%function [LL_I, LL_E, bb, K_II, K_IE, K_EE, m_I, m_E, m] = ...
%    lin_evol_rb_operators(model, [detailed_data])
%
% function computing the time-dependent reduced basis operators and
% vectors. 
%
%Function supports affine decomposition, i.e. different operation modes
% guided by optional field decomp_mode in params. See also the
% contents.txt for general explanation
%
% Required fields of model
% operators_algorithm: name of function for computing the
%                 L_E,L_I,b-operators with arguments (grid, params)
%                 example fv_operators_implicit
%
% 
% Optional fields of model:
%   mu_names : names of fields to be regarded as parameters in vector mu
%   decomp_mode: operation mode of the function
%     -0= 'complete' (default): no parameter dependence or decomposition is 
%                performed. output is as described above.
%     -1= 'components': For each output argument a cell array of output
%                 arguments is returned representing the q-th component
%                 independent of the parameters given in mu_names  
%     -2= 'coefficients': For each output argument a cell array of output
%                 arguments is returned representing the q-th coefficient
%                 dependent of the parameters given in mu_names  
%
% In 'coefficients' mode the detailed data is empty.

% Bernard Haasdonk 23.7.2006

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

% determine affine_decomposition_mode as integer  
decomp_mode = model.decomp_mode;


if decomp_mode == 0 % complete: simple projection on RB
  
  %eval('u0 = ',params.init_value_algorithm,'(grid,params)']);
  [L_I, L_E, b] = model.operators_ptr(model,detailed_data);

  % apply all operator contributions to RB:
  L_E_RB = L_E * detailed_data.RB;
  L_I_RB = L_I * detailed_data.RB;
  
  % fill output quantities
  %    LL_I    :  matrix 
  A = detailed_data.W;
  
  LL_I = detailed_data.RB' * A * L_I_RB;
  
  %    LL_E    :  matrix 
  LL_E = detailed_data.RB' * A * L_E_RB;
  
  %    bb      :  vector 
  bb = detailed_data.RB' * A * b;
  
  %    K_II   : matrix 
  K_II = L_I_RB' * A * L_I_RB;
 
  %    K_IE   : matrix 
  K_IE =  L_I_RB' * A *  L_E_RB;
    
  %    K_EE   :  matrice 
  K_EE =  L_E_RB' * A * L_E_RB;
	  
  %    m_I    :  vector 
  m_I  =  L_I_RB' * A * b;
  
  %    m_E    :  vector 
  m_E  = L_E_RB' * A * b;

  %    m      :  scalar
  m  = 	b' * A * b;
  
elseif decomp_mode == 1
  
  %eval('u0 = ',params.init_value_algorithm,'(grid,params)']);
  [L_I, L_E, b] = model.operators_ptr(model, detailed_data);
  
  Q_L_I = size(L_I,1);
  Q_L_E = size(L_E,1);
  Q_b = length(b);

  Nmax = size(detailed_data.RB,2);
  H = size(detailed_data.RB,1);
  A = detailed_data.W;
  
  % init output quantities
  % scalars:
  m = cell(Q_b*Q_b,1);
  m(:) = {zeros(1,1)};
  % vectors:
  m_I = cell(Q_L_I*Q_b,1);
  m_I(:) = {zeros(Nmax,1)};
  m_E = cell(Q_L_E*Q_b,1);
  m_E(:) = {zeros(Nmax,1)};
  bb = cell(Q_b,1);
  bb(:) = {zeros(Nmax,1)};
  % matrices:
  K_II = cell(Q_L_I*Q_L_I,1);
  K_II(:) = {zeros(Nmax,Nmax)};
  K_IE = cell(Q_L_I*Q_L_E,1);
  K_IE(:) = {zeros(Nmax,Nmax)};
  K_EE = cell(Q_L_E*Q_L_E,1);
  K_EE(:) = {zeros(Nmax,Nmax)};
  LL_I = cell(Q_L_I,1);
  LL_I(:) = {zeros(Nmax,Nmax)};
  LL_E = cell(Q_L_E,1);
  LL_E(:) = {zeros(Nmax,Nmax)};
  
  % apply all operator contributions to RB:
  L_E_RB = cell(Q_L_E,1);
  L_E_RB(:) = {zeros(H,Nmax)};
  L_I_RB = cell(Q_L_I,1);
  L_I_RB(:) = {zeros(H,Nmax)};
  for q=1:Q_L_E
    L_E_RB{q}(:,:)= L_E{q} * detailed_data.RB;
  end;
  for q=1:Q_L_I
    L_I_RB{q}(:,:)= L_I{q} * detailed_data.RB;
  end;  
  
  % fill output quantities
  %    LL_I    :  matrices 
  for q = 1:Q_L_I;
    LL_I{q} = detailed_data.RB' * A * L_I_RB{q};
  end;
  
  %    LL_E    :  matrices 
  for q = 1:Q_L_E;
    LL_E{q} =  detailed_data.RB' * A * L_E_RB{q};
  end;
  
  %    bb      :  vectors 
  for q = 1:Q_b;
    bb{q} = detailed_data.RB' * A * b{q};
  end;

  %    K_II   : matrices 
  for q1 = 1:Q_L_I;
    for q2 = 1:Q_L_I;
      K_II{ (q1-1)*Q_L_I + q2}(:,:) = ... 
	  L_I_RB{q1}' * A * L_I_RB{q2};
    end;
  end;
  
  %    K_IE   : time sequence of matrices 
  % K_IE{1} = Gram matrix of L_I_RB{1} and L_E_RB{1}
  % K_IE{2} = Gram matrix of L_I_RB{1} and L_E_RB{2} 
  % so L_E is "inner loop". This must fit to the sigma-ordering
  for q1 = 1:Q_L_I;
    for q2 = 1:Q_L_E;
      K_IE{ (q1-1)*Q_L_E + q2} = ...
	  L_I_RB{q1}' * A * L_E_RB{q2};
    end;
  end;
  
  %    K_EE   :  matrices 
  for q1 = 1:Q_L_E;
    for q2 = 1:Q_L_E;
      K_EE{ (q1-1)*Q_L_E + q2} = ...
	  L_E_RB{q1}' * A * L_E_RB{q2};
    end;
  end;
  
  %    m_I    :  vectors 
  for q1 = 1:Q_L_I;
    for q2 = 1:Q_b;
      m_I{ (q1-1)*Q_b + q2} = ...
	  L_I_RB{q1}' * A *  b{q2};
    end;
  end;
  
  %    m_E    :  vectors 
  for q1 = 1:Q_L_E;
    for q2 = 1:Q_b;
      m_E{ (q1-1)*Q_b + q2} = ...
	  L_E_RB{q1}' * A * b{q2};
    end;
  end;
  
  %    m      :  scalars 
  for q1 = 1:Q_b;
    for q2 = 1:Q_b;
      m{ (q1-1)*Q_b + q2} = ...
	  b{q1}' * A * b{q2};
    end;
  end;
 
else % decomp_mode== 2 -> coefficients 
  [L_I, L_E, b] = model.operators_ptr(model,[]);
  
  Q_L_I = length(L_I(:));
  Q_L_E = length(L_E(:));
  Q_b = length(b(:));
  
  % rb-operator coefficients can simply be forwarded
  LL_E = L_E;
  LL_I = L_I;
  bb = b;
  
  % pairs: products of sigmas:
  K_II = repmatrows(L_I(:), Q_L_I) .*  repmat(L_I(:), Q_L_I, 1); 
  K_IE = repmatrows(L_I(:), Q_L_E) .*  repmat(L_E(:), Q_L_I, 1); 
  K_EE = repmatrows(L_E(:), Q_L_E) .*  repmat(L_E(:), Q_L_E, 1); 
  m_I  = repmatrows(L_I(:), Q_b)   .*  repmat(b(:)  , Q_L_I, 1); 
  m_E  = repmatrows(L_E(:), Q_b)   .*  repmat(b(:)  , Q_L_E, 1); 
  m    = repmatrows(b(:)  , Q_b)   .*  repmat(b(:)  , Q_b  , 1); 
  
end;

