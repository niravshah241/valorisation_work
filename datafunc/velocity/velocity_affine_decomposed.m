function [vel,lambda] = velocity_affine_decomposed(glob, model)
%function [vel,lambda] = velocity_affine_decomposed(glob, model)
%
% function evaluating a function in the list of global coordinates
% specified in the columns of glob. Result is a matrix of velocity
% vectors as columns of vel.
%
% the function is given by function pointer to components and coefficients:
%    model.velocity_coefficients_ptr
%                  Function must
%                  allow call U = f(model)
%                  Return is a vector of length Q with (parameter
%                  dependent) coefficients in U.
%    model.velocity_components_ptr 
%                  Functions must
%                  allow call U = f(glob,model)
%                  Return of components is a cell array of matrices of
%                  the same size as glob with the point values
%                  of the components. 
%
% Linear combination of components by coefficients then yields the
% complete evaluation.

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

% glob column check
if model.debug
  if ~isempty(glob) && size(glob,1) < size(glob,2)
    warning('coordinates in variable glob are given row-wise, but expected them to be column-wise');
    if model.debug > 2
      keyboard;
    end
  end
end

% determine affine_decomposition_mode as integer
decomp_mode = model.decomp_mode;

if decomp_mode == 2
  vel = model.velocity_coefficients_ptr(model);
  %| \todo lambda needs to be set to something more
  % reasonable for reduced simulations
  lambda = 0;
elseif decomp_mode == 1; % components
  vel = model.velocity_components_ptr(glob,model);
  %| \todo lambda needs to be set to something more
  % reasonable for reduced simulations
  lambda = 0;
else % decomp_mode = 0, complete
  vcoefficients = model.velocity_coefficients_ptr(model);
  vcomponents = model.velocity_components_ptr(glob,model);
  vel = lincomb_sequence(vcomponents,vcoefficients);
  lambda = 1/max(sqrt(vel(:,1).^2 + vel(:,2).^2));
end;

 
