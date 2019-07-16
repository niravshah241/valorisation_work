function [LL_I, bb_I] = nonlin_evol_rb_operators(model, detailed_data, reduced_data)
%function [LL_I, bb_I] = nonlin_evol_rb_operators(model, detailed_data, reduced_data)
%
% function computing the time-dependent reduced basis operator and
% vector. If the decomposition mode is 'coefficients', the detailed
% data are superfluous, can (and should for H-independence) be empty.
%
% required fields of model
% implicit_operators_algorithm: name of function for computing the
%                 L_I,b_I-operators with arguments (grid, params)
%                 example fv_operators_implicit
% inner_product_name: name of inner-product function for
%                     discrete function l2-scalar-product
%
% Function supports affine decomposition, i.e. different operation modes
% guided by optional field affine_decomp_mode in params. See also the
% contents.txt for general explanation
%
% optional fields of params:
%   mu_names : names of fields to be regarded as parameters in vector mu
%   affine_decomp_mode: operation mode of the function
%     'complete' (default): no parameter dependence or decomposition is
%                performed. output is as described above.
%     'components': For each output argument a cell array of output
%                 arguments is returned representing the q-th component
%                 independent of the parameters given in mu_names
%     'coefficients': For each output argument a cell array of output
%                 arguments is returned representing the q-th coefficient
%                 dependent of the parameters given in mu_names
%
% In 'coefficients' mode the detailed data is empty.

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


% Bernard Haasdonk 23.7.2006

% determine affine_decomposition_mode as integer
%decomp_mode = get_affine_decomp_mode(model);

decomp_mode = model.decomp_mode;

if decomp_mode ==0 % complete: simple projection on RB

  A = detailed_data.W;

  RBtranspA = detailed_data.RB' * A;

%  if isequal(params.model_type, 'implicit_nonaffine_linear')
%    L_I = detailed_data.implicit_operator;
%    b_I = detailed_data.implicit_constant;
%  else
  [L_I, b_I] = model.implicit_operators_algorithm(model, detailed_data);
%  end

  % apply all operator contributions to RB:
  L_I_RB = L_I * detailed_data.RB;

  % fill output quantities
  %    LL_I    :  matrix
  LL_I = RBtranspA * L_I_RB;

  %    bb      :  vector
  bb_I = RBtranspA * b_I;

elseif decomp_mode == 1

  if model.newton_solver || model.implicit_nonlinear
    disp('what is the use of the flag implicit_nonlinear ??');
%    [dummy, b_I{1}] = fv_operators_diff_implicit(model, detailed_data, [], []);
    L_I = cell(0,1);
    b_I = cell(0,1);
  else
    [L_I, b_I] = model.implicit_operators_algorithm(model, ...
                                                    detailed_data);
  end

  Q_L_I = size(L_I,1);
  Q_b_I = length(b_I);

  RB = detailed_data.RB;
  Nmax = size(RB,2);
  A = detailed_data.W;

  RBtranspA = RB' * A;

  % init output quantities
  % vectors:
  bb_I = cell(Q_b_I,1);
  bb_I(:) = {zeros(Nmax,1)};
  % matrices:
  LL_I = cell(Q_L_I,1);
  LL_I(:) = {zeros(Nmax,Nmax)};

  % apply all operator contributions to RB:
  L_I_RB = cell(Q_L_I,1);
  L_I_RB(:) = {zeros(size(RB))};
  for q=1:Q_L_I
    L_I_RB{q}(:,:)= L_I{q} * RB;
  end;

  % fill output quantities
  %    LL_I    :  matrices
  for q = 1:Q_L_I;
    LL_I{q} = RBtranspA * L_I_RB{q};
  end;

  %    bb_I      :  vectors
  for q = 1:Q_b_I;
    bb_I{q} = RBtranspA * b_I{q};
  end;

else % decomp_mode== 2 -> coefficients

  [L_I, b_I] = model.implicit_operators_algorithm(model, []);

  % rb-operator coefficients can simply be forwarded
  LL_I = L_I;
  bb_I = b_I;

end;

%| \docupdate 
