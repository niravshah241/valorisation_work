function errs = fv_error(U1,U2,grid,model)
%function errs = fv_error(U1,U2,grid,model);
%
% function computing the error between the two time-dependent
% fv-functions in the columns of U1,U2. Result is a single value or
% sequence of values representing either the l2 differences of the
% temporal sequences or the energy error.
% correct Omega-integrals are computed by respecting the cell-areas
% defined in grid
%
% required fields of model:
%   error_norm: 'l2' or 'energy'
%
% plus additional fields required by fv_l2_error, fv_energy_error.

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

  
% Bernard Haasdonk 30.8.2006
  
if isequal(model.error_norm,'l2')
  errs = fv_l2_error(U1,U2,grid,model);
elseif isequal(model.error_norm, 'energy')
  errs = fv_energy_error(U1,U2,grid,model);
elseif isequal(model.error_norm, 'l2l2')
  errs = fv_l2l2_error(U1,U2,grid,model);
elseif isequal(model.error_norm, 'l1l2')
  errs = fv_l1l2_error(U1,U2,grid,model);
else
  error('error_norm unknown!');
end;
  
