function l2_error = fv_l2_error(U1,U2,grid,model)
%function l2_error = fv_l2_error(U1,U2,grid,params)
%
% function computing the l2-error between the two fv-functions or function
% sequences in U1,U2. Result is a single value or sequence of values 
% corresponding to the column-differences between U1 and U2
% correct Omega-integrals are computed by respecting the cell-areas
% defined in grid,params. Actually, params is currently
% superfluous, but kept for consistency of command line arguments.

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


% Bernard Haasdonk 20.7.2006
  
%  l2_error = sqrt(sum(((U1-U2).^2).*repmat(grid.A(:),1,size(U1,2))));
A = sparse(1:size(U1,1),1:size(U1,1),grid.A(:));
l2_error = sqrt(sum(A * ((U1-U2).^2)));


