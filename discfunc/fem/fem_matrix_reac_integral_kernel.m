function res = fem_matrix_reac_integral_kernel(x,model,df_info,i,j)
% function res = fem_matrix_reac_integral_kernel(x,model,df_info,i,j)
%
% auxiliary function for integral kernel for A_reac
% integral kernel on all elements simultaneously
% here x is in reference triangle!
% f(x) = c hatphi_j hatphi_i 
% multiplication with |det(DF)| is realized in caller after quadrature.
% function can integrate cell-array valued data

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

hat_phi_i = fem_evaluate_basis_function(df_info,x,i);
hat_phi_j = fem_evaluate_basis_function(df_info,x,j);
c = model.reaction(df_info.grid,1:df_info.grid.nelements,x,model); % local model!
if ~iscell(c)
  % the following is vectorial of size nelements!
  res = hat_phi_j * hat_phi_i * c;
else % iscell!!!
  res = cell(1,length(c));
  for q = 1:length(c)
    res{q} = hat_phi_j * hat_phi_i * c{q};
  end;
end;