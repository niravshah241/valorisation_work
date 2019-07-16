function gneu = neumann_values_affine_decomposed(glob, U, normals, params)
%function gneu = neumann_values_affine_decomposed(glob, U, normals, params)
%
% function evaluating a function in the list of global coordinates
% specified in the columns of glob. Result is a row vector of neumann
% boundary values results as columns of gneu.
%
% the function is given by function pointer to components and coefficients: 
%    params.neumann_values_coefficients_ptr 
%                  Function must
%                  allow call gneu = f(params)
%                  Return is a vector of length Q with (parameter
%                  dependent) coefficients in gneu.
%    params.neumann_values_components_ptr 
%                  Functions must
%                  allow call gnue = f(glob,U,normals,params)
%                  Return of components is a cell array of vectors of
%                  the same size as size(glob,1) with the values 
%                  of the neumann boundary function
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


% Bernard Haasdonk 27.8.2009

% determine affine_decomposition_mode as integer  
% glob column check
% normals column check
if params.debug
  if ~isempty(normals) && size(glob,1) < size(glob,2)
    warning('coordinates in variable normals are given row-wise, but expected them to be column-wise');
    if params.debug > 2
      keyboard;
    end
  end
end
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
  gneu = params.neumann_values_coefficients_ptr(params);
elseif decomp_mode == 1; % components
  gneu = params.neumann_values_components_ptr(glob,params);
else % decomp_mode = 0, complete
  gneu_coefficients = params.neumann_values_coefficients_ptr(params);
  gneu_components = params.neumann_values_components_ptr(glob,params); 
  gneu = lincomb_sequence(gneu_components,gneu_coefficients);	  
end;
%| \docupdate 
