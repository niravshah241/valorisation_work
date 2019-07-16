function l2l2_error = fv_l2l2_error(U1,U2,grid,model)
%function l2l2_error = fv_l2l2_error(U1,U2,grid,model)
%
% function computing the l2-error between the two fv-functions or function
% sequences in U1,U2. Result is a single value or sequence of values 
% corresponding to the column-differences between U1 and U2
% correct Omega-integrals are computed by respecting the cell-areas
% defined in grid,model. Params defines the time-spacing for
% correct time integration
%
% note, that the last snapshot does not contribute to the l2 error,
% as FV functions are assumed to be piecewise constant on a forward 
% time-intervall.

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


% Bernard Haasdonk 7.11.2008
  
%  l2_error = sqrt(sum(((U1-U2).^2).*repmat(grid.A(:),1,size(U1,2))));
A = sparse(1:size(U1,1),1:size(U1,1),grid.A(:));
%linf_l2_error_sqr = sqrt(sum(A * ((U1-U2).^2)));
linf_l2_error_sqr = sum(A * ((U1-U2).^2));
l2l2_error = sqrt(cumsum(linf_l2_error_sqr * model.T/model.nt));
l2l2_error = [0,l2l2_error(1:end-1)];

