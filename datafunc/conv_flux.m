function flux = conv_flux(model,U,X,Y)
%function flux = conv_flux(model,[U],[X],[Y])
%
% function computing the convective flux of a convection problem. 
% flux consists of the fields
% Fx, Fy and lambda, representing the x/y-coordinates of the
% velocity in the edge midpoints and a CFL-bound. 
% lambda is only returned in affine_decomp_mode == 'none', as it
% does not make sense otherwise.
% lambda satisfies
%       lambda * sup_u n_jl * f'(u) <= 1
%    e.g. lambda := 1/sup|v(x,y)|  in case of linear flux and u bounded by 1
%    only reasonable in affine_decomp_mode=='none', otherwise not returned.
%
% required fields of model:
% name_flux: 'parabola', 'circular', etc.
% t: real time value in case of time-dependent flux  
%
% if name_flux == 'burgers_parabola' (reasonable domain [0,1]^2)
%    then a parabolic velocity field going to the right is used with
%    maximum velocity model.vmax and angle model.vrot_angle to
%    the x-axis. The computed flux then is v * u^model.flux_pdeg. The input
%    quantity U is assumed to be bounded by linfty-norm 1
%
% if name_flux == 'parabola' (reasonable domain [0,1]^2)
% c : velocity of parabola-profile in model.yrange, constant in x-direction
%
% if name_flux == 'gdl'   (reasonable domain [0,1e-3]x[0,0.25e-3])
% v_in_max : max velocity  
%            a velocity profile is used, which is computed from an
%            elliptic problem: parabolic inflow profile on the left,
%            parabolic outflow on the right, noflow on bottom and middle
%            of top. Dirichlet values in pressure on two top domains
%            yield the velocity. 
% if name_flux == 'gdl2'   (reasonable domain [0,1e-3]x[0,0.25e-3])
% lambda : factor for velocity computation: v = - lambda * grad(p)  
%            a velocity field is used, which is computed from an
%            elliptic problem for the pressure by compute_pressure_gdl2
%            (some gdl-elements in a row, Neumann-0 and dirichlet
%            conditions between 4 and 1 for the pressure, linear
%            decreasing along the channel)
% if name_flux == 'gdl_circ_par_mix' then a mixture of parabolic and
%            circular hand crafted velocity profile is used.
%            does not work good, produces singularities in fv simulation.    
% if name_flux == 'gdl_pgrad_and_parabola' then a mixture of a parabola
%             and a simulated pressure gradient is used. By this,
%             symmetric boundary conditions can be used
%
% Function supports affine decomposition, i.e. different operation modes
% guided by optional field affine_decomp_mode in model. See also the 
% contents.txt for general explanation
%
% optional fields of model:
%   mu_names : names of fields to be regarded as parameters in vector mu
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
% in 'coefficient' mode, the parameters in brackets are empty

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


% Bernard Haasdonk 4.4.2006

flux = [];

% determine affine_decomposition_mode as integer  
%decomp_mode = get_affine_decomp_mode(model);
decomp_mode = model.decomp_mode;
% flag indicating whether the computation respected the decomposition
respected_decomp_mode = 0;

%disp('halt in conv_flux')
%keyboard;

% keep the following up-to-date with conv_flux_linearization!!!
%linear_fluxes = {'gdl2','parabola','gdl_circ_par_mix','gdl',...
%		 'gdl_pgrad_and_parabola'};

%if ismember(model.name_flux, linear_fluxes)
% if fluxname is linear fluxes
if model.flux_linear
  % evaluate linear flux
  lin_flux = conv_flux_linearization(model,U,X,Y);
  if decomp_mode == 0 % i.e. 'none'
    flux.lambda = lin_flux.lambda; 
    if ~isempty(lin_flux.Vx) % otherwise size 0/0 => 0/1 !!!!!
      flux.Fx = lin_flux.Vx(:).* U(:);
      flux.Fy = lin_flux.Vy(:).* U(:);
    else
      flux.Fx = [];
      flux.Fy = [];
    end;
  elseif decomp_mode == 1 % i.e. 'components'
    flux = cell(size(lin_flux));
    Q = length(lin_flux);
    for q = 1:Q
%      f.lambda = lin_flux{q}.lambda; 
      f.Fx = lin_flux{q}.Vx(:).* U(:);
      f.Fy = lin_flux{q}.Vy(:).* U(:);
      flux{q} = f;
    end;
  elseif decomp_mode == 2 % i.e. 'coefficients'
    flux = lin_flux; % copy of coefficient vector from linear flux
  end;
  respected_decomp_mode = 1;
else
  % if no linear flux, then select appropriate nonlinear flux
  switch model.name_flux
   case 'burgers_parabola'
    % i.e. a parabolic velocity field times u^2, which is rotated by
    % an arbitrary angle around the midpoint (0.5,0.5)
    % u is assumed to be bounded by 1, such that lambda can be
    % specified in advance
    
    % NOTE: KEEP THIS DATA IMPLEMENTATION CONSISTENT WITH ITS DERIVATIVE
    % IN conv_flux_linearization  !!!
    P = [X(:)'-0.5; Y(:)'-0.5];
    Rinv = [cos(-model.vrot_angle), - sin(-model.vrot_angle); ...
	    sin(-model.vrot_angle),   cos(-model.vrot_angle) ];
    RP = Rinv*P; % rotated coordinates as columns
    RP(1,:) = RP(1,:) + 0.5;
    RP(2,:) = RP(2,:) + 0.5;
    V = zeros(size(RP));
    V(1,:) = (1- RP(2,:)) .* RP(2,:) * model.vmax * 4;
    V(2,:) = 0;
    RV = Rinv^(-1) * V;
    flux.Fx = RV(1,:)' .* U(:).^model.flux_pdeg;   
    flux.Fy = RV(2,:)' .* U(:).^model.flux_pdeg;

    umax = 1;
    i = find(abs(U)>1.5);
    if ~isempty(i)
      disp(['U not bounded by 1 as assumed by flux, lambda may be ',...
	    'too large.']);
    end;
    flux.lambda = 1/ (model.vmax * 2 * umax^(model.flux_pdeg-1));
    % NOTE: KEEP THIS DATA IMPLEMENTATION CONSISTENT WITH ITS DERIVATIVE
    % IN conv_flux_linearization  !!!
   case 'burgers'
    % i.e. constant velocity v = (model.flux_vx, model.flux_vy) times u^pdeg
    % u is assumed to be bounded by 1, such that lambda can be
    % specified in advance
    Upower =  real(U(:).^model.flux_pdeg);
    flux.Fx = model.flux_vx * ones(size(U(:))) .* Upower;
    flux.Fy = model.flux_vy * ones(size(U(:))) .* Upower;
    
    umax = 1;
    i = find(abs(U)>1.5);
    if ~isempty(i)
      disp(['U not bounded by 1 as assumed by flux, lambda may be ',...
	    'too large.']);
    end;
    vmax = sqrt(model.flux_vy^2 + model.flux_vx^2);
    flux.lambda = 1/ (vmax * 2 * umax^(model.flux_pdeg-1));
   case 'richards'
    flux.Fx = zeros(size(U));
    flux.Fy = zeros(size(U));

    p_mu    = spline_select(model);
    [ breaks, coeffs, pieces, order ] = unmkpp(p_mu);
    p_mu_d  = mkpp(breaks, coeffs(:,1:order-1) .* repmat([order-1:-1:1],pieces,1));


    denom = 1 + ppval(p_mu, X);
    flux.Fx = - 0.001 * ppval(p_mu_d, X) ./ denom .* U';
    flux.Fy = - 0.001 * ppval(p_mu_d, X).^2 .* Y ./ denom.^2 .* U';
    flux.lambda = 1/ (max(max([flux.Fx, flux.Fy])));

   case 'none'
    flux.Fx = zeros(size(U));
    flux.Fy = zeros(size(U));
    % assumption: vmax = 1, umax = 1;
    vmax = 1;
    umax = 1;
    flux.lambda = 1/ (vmax * 2 * umax);
    
   otherwise
    error('flux function unknown');
  end;
  
end;

if decomp_mode>0 & respected_decomp_mode==0
  error('function does not support affine decomposition!');
end;
 
%| \docupdate 
