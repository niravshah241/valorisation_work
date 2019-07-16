function a0 = rb_init_values_default(model,detailed_data)
%function a0 = rb_init_values_default(models,detailed_data)
%
% function computing the reduced basis initial values. If the
% decomposition mode is 'coefficients', the detailed_data are
% superfluous, can (and should for H-independence) be empty.
%
% Required fields of model 
% init_values_algorithm: name of function for computing the
%                 detailed initvalues-DOF with arguments (grid, params)
%                 example: init_values_cog
%
% Function supports affine decomposition, i.e. different operation modes
% guided by optional field decomp_mode in params. See also the
% contents.txt for general explanation
%

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


% Bernard Haasdonk 23.7.2006

% determine affine_decomposition_mode as integer  
decomp_mode = model.decomp_mode;

if decomp_mode == 0 % complete: simple projection on RB
  Nmax = size(detailed_data.RB,2);
  u0 = model.init_values_algorithm(model,detailed_data);
  A = detailed_data.W;
  a0 = u0' * A * detailed_data.RB;
  %a0 = feval(params.inner_product_name, ...
  %      detailed_data.RB,u0,detailed_data.grid,params);
  %  equivalent for fv: a0 = (RB.*repmat(grid.A(:),1,Nmax))' * u0(:);
elseif decomp_mode == 1
  Nmax = size(detailed_data.RB,2);
  u0 = model.init_values_algorithm(model,detailed_data);
  Q_u0 = length(u0);
  a0 = cell(Q_u0,1);
  a0(:) = {zeros(Nmax,1)}; 
  %    a0     : initial data projected on RB set: starting coefficients
  A = detailed_data.W;
  for q=1:Q_u0
    a0{q} = u0{q}' * A *  detailed_data.RB;
  end;
  if model.debug && length(u0) > 1
    if max(abs(detailed_data.RB * a0{1}'-u0{1}))>0.001 
      U0 = [u0{1},u0{2}]; 
      U0_RB = [detailed_data.RB*a0{1}',detailed_data.RB*a0{2}']; 
      plot_params = [];
      plot_params.plot = model.plot;
      plot_params.title = 'U0 Components';
      plot_sequence(U0,detailed_data.grid,plot_params);
      plot_params.title = 'U0_RB';
      plot_sequence(U0_RB,detailed_data.grid,plot_params);
      error('rb initial data badly approximating initial data!!');
    end;
  end;
%  disp('halt in rb_init_values')
%  keyboard;
else % decomp_mode== 2 -> coefficients simply transfered from u0
  u0 = model.init_values_algorithm(model, []);
  a0 = u0;
end;

