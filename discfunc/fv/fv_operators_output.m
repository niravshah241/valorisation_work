function v = fv_operators_output(model,model_data)
%function v = fv_operators_output(model,model_data)
%
% function returning componetns, coefficients, and complete
% operator for a linear output functional on fv discrete functions
% v * U for a dof vector U gives the desired output of the simulation
% assuming, that the output functional is given by 
%
%  `    int_Omega f u dx `
%
% as an integral of a given analytical function multiplied with the
% discrete function.
%
% Function supports affine decomposition, i.e. different operation modes
% guided by optional field decomp_mode in model. See also the
% contents.txt for general explanation

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


% Bernard Haasdonk 13.7.2006
  
% determine affine_decomposition_mode as integer  
decomp_mode = model.decomp_mode;

grid = [];
if ~isempty(model_data)
  grid = model_data.grid;
end

switch decomp_mode
 case 2
   % simply forward coefficients of output_function
   v = model.output_function_ptr([],model);
 case 1
  % simple quadrature quadrature of evaluation of fv function in the midpoints 
  f = model.output_function_ptr([grid.CX(:),grid.CY(:)],model);
  Q_v = length(f);
  v = cell(Q_v,1);
  for q = 1:length(v)
    v{q} = grid.A(:).*(f{q}(:));
  end;
 case 0
  % simple quadrature quadrature of evaluation of fv function in the midpoints 
  f = model.output_function_ptr([grid.CX(:),grid.CY(:)],model);
  % vector that realizes integral over domain of f times U
  % s = \int_Omega f u = sum_e |e| f(c(e)) u(c(e)) = v * U
  % hence v = (|e| f(c(e)))_e;
  v = grid.A(:).*f(:);
end;


