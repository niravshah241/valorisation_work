function K_one = H10_scalar_product_matrix(p,t,c,a,f)
%this routine computes the stiffness matrix for a constant c=1.
%K_one = int_omega (grad phi_j) . grad phi_i dx

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


c_ones=ones(1,length(c));

[K_one,dummy1,dummy2]=assema(p,t,c_ones,a,f);
