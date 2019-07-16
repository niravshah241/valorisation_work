function res = fem_matrix_adv_integral_kernel(x,model,df_info,i,j)
% function res = fem_matrix_adv_integral_kernel(x,model,df_info,i,j)
%
% auxiliary function for integral kernel for A_adv
% integral kernel on all elements simultaneously
% here x is in reference triangle!
% f(x) = (hatphi_j) b^T * (JIT * grad hatphi_i) 
% minus and multiplication with |det(DF)| is realized in caller 
% after quadrature.
% function can integrate cell array valued functions

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
hat_phi_j = fem_evaluate_basis_function(df_info,x,j);
v = model.velocity(df_info.grid,1:df_info.grid.nelements,x,model);%local model!
JIT = df_info.grid.JIT;

if ~iscell(v)
  % the following is vectorial of size nelements!
  res = hat_phi_j * ( v(:,1) .* (JIT(:,1,1) * gradient_hat_phi_i(1) ...
				 + JIT(:,1,2)* gradient_hat_phi_i(2)) ...		  
		      + v(:,2) .* (JIT(:,2,1) * gradient_hat_phi_i(1) ...
				   + JIT(:,2,2)* gradient_hat_phi_i(2)));
else
  res = cell(1,length(v));
  for q = 1:length(res)
    res{q} = hat_phi_j * ( v{q}(:,1) .* (JIT(:,1,1) * gradient_hat_phi_i(1) ...
				 + JIT(:,1,2)* gradient_hat_phi_i(2)) ...		  
		      + v{q}(:,2) .* (JIT(:,2,1) * gradient_hat_phi_i(1) ...
				   + JIT(:,2,2)* gradient_hat_phi_i(2)));
    
  end;
end;