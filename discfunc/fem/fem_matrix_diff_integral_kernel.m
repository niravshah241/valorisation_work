function res = fem_matrix_diff_integral_kernel(x,model,df_info,i,j)
%function res = fem_matrix_diff_integral_kernel(x,model,df_info,i,j)
%
% auxiliary function for integral kernel for A_diff
% integral kernel on all elements simultaneously
% here x is in reference triangle!
% ``f(x) = (\nabla \hat\phi_j)^T * JIT^T * A^T * (JIT * \nabla \hat\phi_i)``
% multiplication with `|det(DF)|` is realized in caller after quadrature.
% function can integrate cell-array valued function

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


% B. Haasdonk 22.2.2011
% I. Maier 24.03.2011

gradient_hat_phi_i = fem_evaluate_scalar_basis_function_derivative(df_info,x,i)';
gradient_hat_phi_j = fem_evaluate_scalar_basis_function_derivative(df_info,x,j)';
a = model.diffusivity_tensor(df_info.grid,1:df_info.grid.nelements,x,model); % local model!
JIT = df_info.grid.JIT;

if ~iscell(a)
  % the following are vectorial of size nelements!
  % AT = Atransposed = (a1,a2; a3,a4)
  ATJ_11 =  a(:,1) .* JIT(:,1,1) + a(:,2) .* JIT(:,2,1);
  ATJ_21 =  a(:,3) .* JIT(:,1,1) + a(:,4) .* JIT(:,2,1);
  ATJ_12 =  a(:,1) .* JIT(:,1,2) + a(:,2) .* JIT(:,2,2);
  ATJ_22 =  a(:,3) .* JIT(:,1,2) + a(:,4) .* JIT(:,2,2);
  %JTATJ_11 = ...
  %    ( a(:,1) .* JIT(:,1,1) + a(:,3) .* JIT(:,2,1)) .* JIT(:,1,1) + ... 
  %    ( a(:,1) .* JIT(:,1,2) + a(:,3) .* JIT(:,2,2)) .* JIT(:,2,1);
  %JTATJ_12 = ...
  %    ( a(:,1) .* JIT(:,1,1) + a(:,3) .* JIT(:,2,1)) .* JIT(:,1,2) + ... 
  %    ( a(:,1) .* JIT(:,1,2) + a(:,3) .* JIT(:,2,2)) .* JIT(:,2,2);
  %JTATJ_21 = ...
  %    ( a(:,2) .* JIT(:,1,1) + a(:,4) .* JIT(:,2,1)) .* JIT(:,1,1) + ... 
  %    ( a(:,2) .* JIT(:,1,2) + a(:,4) .* JIT(:,2,2)) .* JIT(:,2,1);
  %JTATJ_22 = ...
  %%    ( a(:,2) .* JIT(:,1,1) + a(:,4) .* JIT(:,2,1)) .* JIT(:,1,2) + ... 
  %    ( a(:,2) .* JIT(:,1,2) + a(:,4) .* JIT(:,2,2)) .* JIT(:,2,2);
  JTATJ_11 = JIT(:,1,1).* ATJ_11 + JIT(:,2,1).* ATJ_21;
  JTATJ_12 = JIT(:,1,1).* ATJ_12 + JIT(:,2,1).* ATJ_22;
  JTATJ_21 = JIT(:,1,2).* ATJ_11 + JIT(:,2,2).* ATJ_21;
  JTATJ_22 = JIT(:,1,2).* ATJ_12 + JIT(:,2,2).* ATJ_22;
  res = ...
      JTATJ_11 * gradient_hat_phi_j(1)* gradient_hat_phi_i(1) + ...
      JTATJ_12 * gradient_hat_phi_j(1)* gradient_hat_phi_i(2)+ ...
      JTATJ_21 * gradient_hat_phi_j(2)* gradient_hat_phi_i(1)+ ...
      JTATJ_22 * gradient_hat_phi_j(2)* gradient_hat_phi_i(2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else % iscell! 
  res = cell(1,length(a));
  for q= 1:length(a);
    % the following are vectorial of size nelements!
    % AT = Atransposed = (a1,a2; a3,a4)
    ATJ_11 =  a{q}(:,1) .* JIT(:,1,1) + a{q}(:,2) .* JIT(:,2,1);
    ATJ_21 =  a{q}(:,3) .* JIT(:,1,1) + a{q}(:,4) .* JIT(:,2,1);
    ATJ_12 =  a{q}(:,1) .* JIT(:,1,2) + a{q}(:,2) .* JIT(:,2,2);
    ATJ_22 =  a{q}(:,3) .* JIT(:,1,2) + a{q}(:,4) .* JIT(:,2,2);
    %JTATJ_11 = ...
    %    ( a(:,1) .* JIT(:,1,1) + a(:,3) .* JIT(:,2,1)) .* JIT(:,1,1) + ... 
    %    ( a(:,1) .* JIT(:,1,2) + a(:,3) .* JIT(:,2,2)) .* JIT(:,2,1);
    %JTATJ_12 = ...
    %    ( a(:,1) .* JIT(:,1,1) + a(:,3) .* JIT(:,2,1)) .* JIT(:,1,2) + ... 
    %    ( a(:,1) .* JIT(:,1,2) + a(:,3) .* JIT(:,2,2)) .* JIT(:,2,2);
    %JTATJ_21 = ...
    %    ( a(:,2) .* JIT(:,1,1) + a(:,4) .* JIT(:,2,1)) .* JIT(:,1,1) + ... 
    %    ( a(:,2) .* JIT(:,1,2) + a(:,4) .* JIT(:,2,2)) .* JIT(:,2,1);
    %JTATJ_22 = ...
    %%    ( a(:,2) .* JIT(:,1,1) + a(:,4) .* JIT(:,2,1)) .* JIT(:,1,2) + ... 
    %    ( a(:,2) .* JIT(:,1,2) + a(:,4) .* JIT(:,2,2)) .* JIT(:,2,2);
    JTATJ_11 = JIT(:,1,1).* ATJ_11 + JIT(:,2,1).* ATJ_21;
    JTATJ_12 = JIT(:,1,1).* ATJ_12 + JIT(:,2,1).* ATJ_22;
    JTATJ_21 = JIT(:,1,2).* ATJ_11 + JIT(:,2,2).* ATJ_21;
    JTATJ_22 = JIT(:,1,2).* ATJ_12 + JIT(:,2,2).* ATJ_22;
    res{q} = ...
	JTATJ_11 * gradient_hat_phi_j(1)* gradient_hat_phi_i(1) + ...
	JTATJ_12 * gradient_hat_phi_j(1)* gradient_hat_phi_i(2)+ ...
	JTATJ_21 * gradient_hat_phi_j(2)* gradient_hat_phi_i(1)+ ...
	JTATJ_22 * gradient_hat_phi_j(2)* gradient_hat_phi_i(2);
  end;
end;

  