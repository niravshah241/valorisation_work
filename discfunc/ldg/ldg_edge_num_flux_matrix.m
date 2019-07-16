function edge_num_flux_mat = ldg_edge_num_flux_matrix(edge_flux_matrix,params)
%function edge_num_flux_mat = ldg_edge_num_flux_matrix(edge_flux_matrix,params)
%
% function computing an edge_num_flux_matrix (or its components)
% from the analytical edge flux matrix edge_flux_matrix, which is
% computed in ldg_operators_adv_explicit
%
% currently simple upwind flux is computed
%
% In "components" mode:
%    input  edge_flux_matrix(i,j,e,f,q): edge integral contribution 
%   and
%    output edge_num_flux_matrix(...): ...
%
% In "complete" mode:
%    input  edge_flux_matrix(i,j,e,f): edge integral contribution 
%   and
%    output edge_num_flux_matrix(...): ...

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

disp(' to be adjusted!!  ');
keyboard;

decomp_mode = params.decomp_mode;

if decomp_mode == 2
  % for upwind: simple coefficients of velocity.
  edge_num_flux_mat = params.velocity_coefficients_fct(params);
elseif decomp_mode ==1
  error('to be implemented!');
  % identify consistent upwind direction from components 

  % sign of diagonal entries are sign of v*n
  
  % find consistent upwind directions
  % => check of affine parameter dependence
  ...
  
  % remove dirichlet inflow contributions

  % replace downwind flux values with upwind flux values

  % return matrix
  
else % "complete"
  disp('complete evaluation should be accelerated in ldg_edge_num_flux_matrix');
  params.decomp_mode = 2;
  edge_num_flux_mat_coeff = ldg_edge_num_flux_matrix(edge_flux_matrix,params);
  params.decomp_mode = 1;
  edge_num_flux_mat_comp = ldg_edge_num_flux_matrix(edge_flux_matrix,params);
  params.decomp_mode = 0;
  edge_num_flux_mat = lincomb_sequence(edge_num_flux_mat_comp,...
				       edge_num_flux_mat_coeff);
end;
%| \docupdate 
