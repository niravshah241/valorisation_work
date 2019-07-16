function res = disc_init_values(model,model_data)
%function res = disc_init_values(model,model_data)
%
% function computing the init-value discrete function by
% l2projection of the init value function on the discrete function space
%
% in 'complete' mode a singe discrete function is produced,
% in 'components' mode, it is a cell array of discrete functions
% in 'coefficients' mode it is a vector of coefficients.

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


% Bernard Haasdonk 22.7.2006

if nargin ~= 2
  error('wrong number of parameters!');
end;

% affine decomposition is supported by init_data, so automatically
% satisfied.

decomp_mode = model.decomp_mode;

if decomp_mode == 2 % = coefficients
  res = model.init_values_ptr([],model);
elseif decomp_mode == 1 % = components
  % idea for improvement: method l2_project sequence: integration of a sequence
  % of functions  
  
  % the following is an expensive workaround
  params = model; 
  params.nelements = model_data.grid.nelements; 
  params.pdeg = model.pdeg;
  % determine dimrange
  tmp = model.init_values_ptr([0,0],model);
  params.dimrange = size(tmp{1},1); % possibly vectorial function

  Q_u0 = length(tmp);
  res = cell(1,Q_u0);
  for q = 1:Q_u0
    params.single_comp_index = q;
    f = @(einds,loc,grid,params) get_single_comp(...
	model.init_values_ptr(...
	    local2global(grid,einds,loc,params),...
	    params),...
	params); 
    res{q} = model.l2project(f,model.init_values_qdeg,model_data.grid, ...
			     params);
  end;
  
else % decomp_mode = 0 == complete
  params = model; 
  params.nelements = model_data.grid.nelements; 
  params.pdeg = model.pdeg;
  % determine dimrange
  tmp = model.init_values_ptr([0,0],model);
  params.dimrange = size(tmp,1); % possibly vectorial function
  f = @(einds,loc,grid,params) model.init_values_ptr(...
      local2global(grid,einds,loc,params),params); 
  res = model.l2project(f,model.init_values_qdeg,model_data.grid, ...
			params);
end;

% extract single array from cell array
function res = get_single_comp(cellarr,params)
res = cellarr{params.single_comp_index};
%| \docupdate 
