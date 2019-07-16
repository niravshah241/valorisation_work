function U0 = init_values_affine_decomposed(glob, params)
%function U0 = init_values_affine_decomposed(glob, params)
%
% function evaluating a function in the list of global coordinates
% specified in the columns of glob. Result is a matrix of initial
% data results as columns of U0.
%
% the function is given by function pointer to components and coefficients: 
%    params.init_values_coefficients_ptr 
%                  Function must
%                  allow call U = f(params)
%                  Return is a vector of length Q with (parameter
%                  dependent) coefficients in U.
%    params.init_values_components_ptr 
%                  Functions must
%                  allow call U = f(glob,params)
%                  Return of components is a cell array of matrices of
%                  the same size as glob with the point values
%                  of the components. 
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


% Bernard Haasdonk 27.8.2009

% determine affine_decomposition_mode as integer  
% glob column check
if params.debug
  if ~isempty(glob) && size(glob,1) < size(glob,2)
    warning('coordinates in variable glob are given row-wise, but expected them to be column-wise');
    if params.debug > 2
      keyboard;
    end
  end
end
decomp_mode = params.decomp_mode;

if decomp_mode == 2
  U0 = params.init_values_coefficients_ptr(params);
elseif decomp_mode == 1; % components
  U0 = params.init_values_components_ptr(glob,params);
else % decomp_mode = 0, complete
  Ucoefficients = params.init_values_coefficients_ptr(params);
  Ucomponents = params.init_values_components_ptr(glob,params); 
  U0 = lincomb_sequence(Ucomponents,Ucoefficients);	  
end;    
%| \docupdate 
