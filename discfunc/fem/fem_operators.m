function [A,r] = fem_operators(model,model_data)
%function [A,r] = fem_operators(model,model_data)
%
% function computing the matrix and rhs of an elliptic problem with
% finite element discretization
% supports affine decomposition, i.e. result 
% depending on model.decomp_mode

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

if model.decomp_mode == 2
  df_info = []; % not required
else
  df_info = model_data.df_info;  
end;

[r_source, r_dirichlet, r_neumann, r_robin] = ...
    fem_rhs_parts_assembly(model,df_info);

[A_diff , A_adv, A_reac, A_dirichlet, A_neumann, A_robin] = ...
    fem_matrix_parts_assembly(model,df_info);  


if model.decomp_mode == 0 % == complete: simple addition

  % assemble right hand side
  r = r_source + r_neumann + r_robin + r_dirichlet;
  if model.verbose > 5
    disp('rhs assembled');
  end;
  % sparse system matrix:
  A = spalloc(model_data.df_info.ndofs,model_data.df_info.ndofs,10); % rough upper bound for number of nonzeros
  A = A_diff + A_adv + A_reac + A_neumann + A_robin + A_dirichlet;

  if model.verbose > 5
    disp('matrix assembled');
  end;
  
elseif model.decomp_mode == 1 % == components: merge to cell arrays

  r = [r_source(:)', r_dirichlet(:)', r_neumann(:)', r_robin(:)'];

  A = [A_diff(:)' , A_adv(:)', A_reac(:)', ...
       A_dirichlet(:)', A_neumann(:)', A_robin(:)'];
  
%  keyboard;
  
else % decomp_mode == 2, coefficients: merge coeff vectors
  
  r = [r_source(:); r_dirichlet(:); r_neumann(:); r_robin(:)];
  
  A = [A_diff(:) ; A_adv(:); A_reac(:); ...
       A_dirichlet(:); A_neumann(:); A_robin(:)];
  
%  keyboard;
  
end;