function diffusivity = diffusivity(model,X,Y,U)
%function diffusivity = diffusivity(model,[X],[Y],[U])
%
% TODO: Add the function values U to parameter list
% function computing the diffusivity of a convection problem by pointwise
% evaluation in the point sequences indicated by X and Y. 
% fields of diffusivity:
%     epsilon: upper bound on diffusivity value
%     K: vector with diffusivity values
%
% fields of model:
% name_diffusivity = 'none', 'homogeneous', 'linear_gradient'
%
% in case of 'homogeneous'
% k : scalar diffusivity parameter
%
% in case of 'linear_gradient'
% diff_left : value at left border
% diff_right : value at right border
%
% optional: time
%
% Function supports affine decomposition, i.e. different operation modes
% guided by optional field affine_decomp_mode in model. See also the 
% contents.txt for general explanation
%
% optional fields of model:
%   affine_decomp_mode: operation mode of the function
%     'none' (default): no parameter dependence or decomposition is 
%                performed. output is as described above.
%     'components': For each output argument a cell array of output
%                 arguments is returned representing the q-th component
%                 independent of the parameters given in mu_names  
%     'coefficients': For each output argument a cell array of output
%                 arguments is returned representing the q-th coefficient
%                 dependent of the parameters given in mu_names  
%
% in 'coefficients' mode, the parameters in brackets are empty

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


% Bernard Haasdonk 11.4.2006

% determine affine_decomposition_mode as integer  
%decomp_mode = get_affine_decomp_mode(model);
decomp_mode = model.decomp_mode;
% flag indicating whether the computation respected the decomposition
respected_decomp_mode = 0;

diffusivity = [];
diffusivity.epsilon = 0;

if isequal(model.name_diffusivity,'none')
    diffusivity.K = zeros(length(X),1);
elseif isequal(model.name_diffusivity,'homogeneous')
  if decomp_mode == 0
    diffusivity.epsilon = model.k;
    diffusivity.K = model.k * ones(length(X),1);
  elseif (decomp_mode == 1)
    d = [];
    % single component independent whether k in mu
    d.epsilon = model.k;
    d.K = ones(length(X),1);
    diffusivity = {d};
  else % decomp_mode ==2, single component
    diffusivity = model.k;
  end;
  respected_decomp_mode = 1;
elseif isequal(model.name_diffusivity,'homogeneous2') ...
        || isequal(model.name_diffusivity,'linear_gradient')

  if decomp_mode == 0
    diffusivity.epsilon = model.k;
    diffusivity.K = model.k * speye(2*length(X));
  elseif (decomp_mode == 1)
    d = [];
    % single component independent whether k in mu
    d.epsilon = model.k;
    d.K = speye(2*length(X));%ones(length(X),1);
    diffusivity = {d};
  else % decomp_mode ==2, single component
    diffusivity = model.k;
  end;
  respected_decomp_mode = 1;
elseif isequal(model.name_diffusivity,'linear_gradient2')
  temp_m = (model.diff_right-model.diff_left) / ...
           (model.xrange(2) - model.xrange(1));
  lin_grad = @(x) (x - model.xrange(1) );
  if decomp_mode == 0
    diffusivity.epsilon = max(model.diff_left, model.diff_right);
    diffusivity.K = model.diff_left + (lin_grad(X) * temp_m);
  elseif (decomp_mode == 1)
    d = [];
    d.epsilon = max(model.diff_left, model.diff_right);
    d.K = lin_grad(X);
    diffusivity = {d};
  else % decomp_mode == 2
    diffusivity = temp_m;
  end
elseif isequal(model.name_diffusivity,'richards_nonlinear')
  % How can I construct a tensor for each dof?
  vlen = size(U,1);
  p_mu = spline_select(model);
  U = U - model.gravity * Y.*(1+ppval(p_mu, X));
  Kaval = model.k * model.richards_k(U') .* model.richards_p(U');
  diffusivity.K = spdiags(reshape([Kaval;Kaval],2*vlen,1),0,2*vlen,2*vlen);
  diffusivity.epsilon = max(Kaval);
  %    Utemp = U' * model.k + 0.0004;
  %    diffusivity.K = spdiags(reshape([Utemp;zeros(1,vlen)],2*vlen,1),0,2*vlen,2*vlen);
  %    diffusivity.epsilon = max(Utemp);
else 
  error('diffusivity function unknown');
end;

if decomp_mode>0 && respected_decomp_mode==0
  error('function does not support affine decomposition!');
end;

%| \docupdate 
