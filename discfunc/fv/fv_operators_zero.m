function [L,b] = fv_operators_zero(model,model_data, NU_ind)
%function [L,b] = fv_operators_zero(model,model_data[, NU_ind])
%
% function returning zero componetns, coefficients, and complete
% zero matrix for use if neumann/diff/conv-imlicit/explicit is to
% be used as zero.

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


% Function supports affine decomposition, i.e. different operation modes
% guided by optional field decomp_mode in model. See also the
% contents.txt for general explanation

% Bernard Haasdonk 13.7.2006
  
% determine affine_decomposition_mode as integer  
decomp_mode = model.decomp_mode;

if nargin < 3
  NU_ind = [];
end

grid = [];
if ~isempty(model_data)
  grid = model_data.grid;
end

switch decomp_mode
 case 2
  L = [];
  b = [];
 case 1
  L = {};
  b = {};
 case 0
  n = grid.nelements;
  if ~isempty(NU_ind)
    m = length(NU_ind);
  else
    m = n;
  end
  L = sparse(m,n);
  b = zeros(m,1);
end;

