function res = fem_h1_norm(df)
%function res = fem_h1_norm(df)
%
% function computing the H1-norm of a discrete fem function

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


ressqr = df.dofs(:)' * ...
	 df.df_info.h10_inner_product_matrix * ...
	 df.dofs(:)+
         df.dofs(:)' * ...
	 df.df_info.l2_inner_product_matrix * ...
	 df.dofs(:);
if ressqr < 0
  ressqr = 0;
end;
res = sqrt(ressqr);
