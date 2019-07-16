function flux_lin = velocity_to_be_distributed_in_single_files(model,U,X,Y)
%function flux_lin = conv_flux_linearization(model,[U],[X],[Y])
%
% function computing the derivative of a convective flux by forward 
% difference or exact evaluation. flux_lin consists of the fields
% Vx, Vy, representing the x/y-coordinates of the
% velocity in the edge midpoints. 
% lambda : bound such that
%      lambda * sup_u n_jl * f'(u) <= 1
%    e.g. lambda := 1/sup|v(x,y)|
%    only reasonable in affine_decomp_mode=='none', otherwise not returned.
%
% For a linear flux in conv_flux the following is identical
%   
%     (Fx(U),Fy(U)) = (Vx,Vy) * U
%
% If the flux is linear, exact evaluation is performed, U is not used
% If the field use_velocity_matrixfile is set, then a filename in 
% 'velocity_matrixfile' is expected, where the point evaluation of
% the velocity are assumed to be stored. 
% The file is expected to be located in [rbmatlabhome,'datafuncs/mat']. 
% The file is expected to have
% variable-vectors Vx, Vy, X, Y, and a scalar lambda. The name_flux
% parameter is in this case ignored (might be used during generation of
% the matrixfile by gen_velocity_matrixfile)  
% As access to this matrix is always depending on
% expensive point-correspondence searches, a file-caching of the 
% velocity extraction can be turned on by setting the field
% filecache_velocity_matrixfile_extract = 1
% This is useful, if the same point-lists X and Y are requested
% multiply, e.g. in timedependent case with constant velocity.
% Still this multiple filecaching can be expensive, instead by
% filecache_velocity_matrixfile_extract = 2, it is assumed, that
% the given file can directly be taken as velocities, i.e. the
% points coincide. For this a selection of the filename has to be
% done by the calling numerical routines.
%
% required fields of model as desired by conv_flux
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
% In 'coefficient' mode, the parameters in brackets are empty

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


% Bernard Haasdonk 21.4.2006
  
flux_lin = [];

% determine affine_decomposition_mode as integer  
%decomp_mode = get_affine_decomp_mode(model);
decomp_mode = model.decomp_mode;
% flag indicating whether the computation respected the decomposition
respected_decomp_mode = 0;

% keep the following up-to-date with conv_flux!!!
%linear_fluxes = {'gdl2','parabola','gdl_circ_par_mix','gdl',...
%		 'gdl_pgrad_and_parabola'};

%if ismember(model.name_flux, linear_fluxes)
if model.flux_linear
  % check if flux is to be computed or read from file
  if isfield(model,'use_velocity_matrixfile') & ...
	(model.use_velocity_matrixfile==1)
    
    % affine decomposition is trivial: no dependency on mu
    
    if decomp_mode < 2 % i.e. 'none' or 'components'
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% loading flux matrix.
      % load V.X: list of X-Koordinaten
      % load V.Y: list of Y-Koordinaten
      % load V.Vy: list of Y-component of velocity 
      % load V.Vx: list of Y-component of velocity 
      % load V.lambda: lambda for CFL computation 
      
      if ~isfield(model,'velocity_matrixfile')
	error('field velocity_matrixfile not set!!');
      end;
      
      fullfn = fullfile(rbmatlabhome,'datafunc','mat',...
			model.velocity_matrixfile); 
      
      if ~exist(fullfn,'file') & ~cache('exist',fullfn)
	error(['warning: velocity_matrixfile not existing. ',...
	       'Call gen_velocity_matrixfile with suitable' ...
	       ' divclean_mode.']);
      end;
       
      % the following file extraction can be expensive, so the
      % extraction can additionally be cached.
      if ~isfield(model,'filecache_velocity_matrixfile_extract')
	model.filecache_velocity_matrixfile_extract = 0;
      end;

      if ~model.filecache_velocity_matrixfile_extract % expensive call
	[flux_lin.Vx, flux_lin.Vy, lambda] = ...
	    velocity_matrixfile_extract(fullfn,X,Y);
      elseif model.filecache_velocity_matrixfile_extract==1 
	% cached functioncall
	[flux_lin.Vx, flux_lin.Vy, lambda] = ...
	    filecache_function('velocity_matrixfile_extract', ...
			       fullfn,X,Y);
      else % mode 2, i.e assume, that velocityfilename is set correctly to
           % requested points
	   %	   V = load_cached(fullfn);
	   V = cache('load',fullfn);
	   % find integer indices such that V.X(j)==X(i) and V.Y(j) = Y(i) 
	   % if whole flux-matrix is requested
	   X = X(:); Y = Y(:);
	   if ~isequal(X,V.X(:)) | ~isequal(Y,V.Y(:))
	     error('matrixfile does not fit to requested points!');
	   end;
	   i = 1:length(X);
	   flux_lin.Vx = V.Vx;
	   flux_lin.Vy = V.Vy;
	   %NaN*ones(size(X));
%	   flux_lin.Vy = NaN*ones(size(X));
%	   flux_lin.Vx(i) = V.Vx(j);
%	   flux_lin.Vy(i) = V.Vy(j);
	   lambda = V.lambda;
      end;
      
      if decomp_mode == 0
	flux_lin.lambda = lambda;
      end;
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end of loading flux matrix.
      if decomp_mode == 1 % put single matrix into cell array
	flux_lin = {flux_lin};
      end;
    else % decomp_mode = 2;
      flux_lin = 1; % factor one for single component
    end;
    respected_decomp_mode = 1;
    
  else %%%%%%%%%%% compute the fluxes instead of loading matrix from file  
    
    if isequal(model.name_flux,'parabola')
      % determine lambda_jl to be globally constant such that
      %   lambda_jl * sup_u n_jl * f'(u) <= 1
      %  f'(u) = v(x,y)
      %  e.g. lambda := 1/sup|v(x,y)|
      if decomp_mode == 0
	flux_lin.lambda = 1/(abs(model.c)/4.0+ 1e-10); % some eps for divby0
      end;
      
      % flux.precomputable = 1; % velocity field can be computed in advance 
      % flux values : velocity field in the edge midpoints
      if decomp_mode==2
	flux_lin = model.c;
      elseif decomp_mode == 1
        dummy.Vy = zeros(length(X),1);
	dummy.Vx = Y(:) .* (1-Y(:));
	flux_lin = {dummy};
      elseif decomp_mode == 0 % decomp_mode 0
	flux_lin.Vy = zeros(length(X),1);
	flux_lin.Vx = model.c * Y(:) .* (1-Y(:));  
      else
	error('unknown decmposition mode!!');
      end;
      respected_decomp_mode = 1;
    elseif isequal(model.name_flux,'gdl2')
      % load pressure from file and generate velocity field
      p = load(fullfile(rbmatlabhome,'datafunc','mat','gdl_pressure2.mat'),...
	       'p','model');
      [ux, uy] = evaluate_gradient(X, Y,p.p,p.model);  
      flux_lin.Vx = -ux * model.lambda;
      flux_lin.Vy = -uy * model.lambda;
      % inverse of maximum absolute velocity value 
      % correct up to a constant factor
      flux_lin.lambda = max([flux_lin.Vx(:); flux_lin.Vy(:)])^-1; 
      %    flux.Fx = ux.* U(:);
      %    flux.Fy = uy.* U(:);
    elseif isequal(model.name_flux,'gdl_circ_par_mix')
      % determine lambda_jl to be globally constant such that
      %   lambda_jl * sup_u n_jl * f'(u) <= 1
      %  f'(u) = v(x,y)
      %  e.g. lambda := 1/sup|v(x,y)|
      
      flux_lin.lambda = 1/abs(model.v_max); % some eps for preventing divby0 
      
      % flux values : parabolic and circular profile
      Vx1 = model.v_max * (100e-6)^(-2) * Y(:) .* (200e-6 -Y(:)); 
      Vy1 = zeros(length(X),1);
      
      % circular field by squared affine radial distance to M = (500,500)*1e-6
      e = (X(:)-500e-6).^2 + (2.5*Y(:)-500e-6).^2; 
      % desired velocity profile: 
      %  vabs(e) = linearly growing from 0 to vmax for 250<=e<=500 
      %  vabs(e) = linearly decreasing to 0 for 500<=e<=750 
      %  vabs(e) = zero otherwise 
      vabs = zeros(size(e));
      se = sqrt(e);
      i = find( se>= 250e-6 & se <= 500e-6);
      vabs(i) = model.v_max / 250e-6 * (se(i)-250e-6);
      i = find( se>= 500e-6 & se <= 750e-6);
      vabs(i) = -model.v_max / 250e-6 * (se(i)-750e-6);
      
      % compute velocity directions by gradient of elliptic field
      dx_e = 2*(X(:)-500e-6);
      dy_e = 5*(Y(:)-500e-6);
      norm_grad_e_inv = (dx_e.^2 + dy_e.^2).^(-1);
      % rotate by 90ø and normalization 
      Vx2 = - dy_e .* norm_grad_e_inv;  
      Vy2 = dx_e .*norm_grad_e_inv;
      % scale to desired absolute velocity
      Vx2 = Vx2 .* vabs;
      Vy2 = Vy2 .* vabs;
      
      flux_lin.Vx = (model.v_weight * Vx1 + (1-model.v_weight) * Vx2);
      flux_lin.Vy = (model.v_weight * Vy1 + (1-model.v_weight) * Vy2);
    elseif isequal(model.name_flux,'gdl')
      % load pressure from file and generate velocity field
      p0 = load('gdl_pressure_v0000.mat','p','model');
      p1000 = load('gdl_pressure_v1000.mat','p','model');
      factor = 1000; % scaling: maximum inflow now 1000/factor m/s
      [ux0, uy0] = evaluate_gradient(po0.model,X, Y,p0.p);  
      ux0 = -ux0/factor;
      uy0 = -uy0/factor;
      [ux1000, uy1000] = evaluate_gradient(p1000.model,X, Y,p1000.p);  
      ux1000 = -ux1000/factor;
      uy1000 = -uy1000/factor;
      % linear interpolation between velocity profiles yields exact solution
      % of corresponding elliptic boundary value problem due to 
      % linearity of the problem.
      ux = model.v_in_max * ux1000 + (1-model.v_in_max) * ux0;
      uy = model.v_in_max * uy1000 + (1-model.v_in_max) * uy0;
      
      % inverse of maximum absolute velocity value 
      flux_lin.lambda = max([ux(:); uy(:)])^-1; 
      % correct up to a constant factor
      flux_lin.Vx = ux;
      flux_lin.Vy = uy;
    elseif isequal(model.name_flux,'gdl_pgrad_and_parabola')
      % load pressure from file and generate velocity field
      p0 = load('gdl_pressure_dirichlet.mat','p','model');
      factor = 1000; % scaling: maximum inflow now 1000/factor m/s
      [ux0, uy0] = evaluate_gradient(p0.model,X, Y,p0.p);  
      ux0 = -ux0/factor;
      uy0 = -uy0/factor;
      uypar = zeros(length(X),1);
      uxpar = model.v_in_max * (model.yrange(2)-model.yrange(1))^-2 * 4 ...
	      * (Y(:)-model.yrange(1)) .* (model.yrange(2)-Y(:));
      
      % linear interpolation between velocity profiles yields exact solution
      % of corresponding elliptic boundary value problem due to 
      % linearity of the problem.
      ux = model.v_in_max * uxpar + (1-model.v_in_max) * ux0;
      uy = model.v_in_max * uypar + (1-model.v_in_max) * uy0;
      
      % inverse of maximum absolute velocity value 
      flux_lin.lambda = max([ux(:); uy(:)])^-1; % correct up to a factor   
      flux_lin.Vx = ux;
      flux_lin.Vy = uy;      
    elseif isequal(model.name_flux, 'transport')
      if decomp_mode == 2
	flux_lin = 1; % single component factor 1
      elseif decomp_mode == 1 % components:
	tmp.Vx = model.transport_x * ones(size(U));
	tmp.Vy = model.transport_y * ones(size(U));
	flux_lin = {tmp};
      else % decomp_mode = 0, complete:
	flux_lin.Vx = model.transport_x * ones(size(U));
	flux_lin.Vy = model.transport_y * ones(size(U));
	flux_lin.lambda = max([flux_lin.Vx(:); 
          		       flux_lin.Vy(:)])^-1; % correct up to a factor
	
      end;
      respected_decomp_mode = 1;
    elseif isequal(model.name_flux, 'richards')

      flux_lin.Vx = zeros(size(U));
      flux_lin.Vy = zeros(size(U));

%      p_mu    = spline_select(model);
%      [ breaks, coeffs, pieces, order ] = unmkpp(p_mu);
%      p_mu_d  = mkpp(breaks, coeffs(1:order-1) .* [order-1:-1:1]);

%      denom = 1 + ppval(p_mu, X);
      [res1, res2] = inv_geo_trans_derivative(model,X,Y,{[1],[2],[1,1],[1,2]},{[1],[2],[2,1],[2,2]});
      d1  = res1{3} + res2{3};
      d2  = res1{4} + res2{4};
      flux_lin.Vx = -model.k .* (res1{1} .* d1 + res1{2} .* d2);
      flux_lin.Vy = -model.k * (res2{1} .* d1 + res2{2} .* d2);
%      flux_lin.Vx = model.k * ppval(p_mu_d, X) ./ denom;
%      flux_lin.Vy = -model.k * ppval(p_mu_d, X).^2 .* Y ./ denom.^2;
      % inverse of maximum absolute velocity value 
      flux_lin.lambda = max([flux_lin.Vx(:); flux_lin.Vy(:)])^-1; % correct up to a factor   
    else
      error('unknown flux!!');
    end; %%% end of flux select
  end;

  
  
  % if not is linear flux: exact computation or forward difference 
elseif isequal(model.name_flux,'burgers_parabola')
  % NOTE: KEEP THIS DATA IMPLEMENTATION CONSISTENT WITH ITS FLUX 
  % IN conv_flux  !!!  
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
  flux_lin.Vx = RV(1,:)' .* model.flux_pdeg .* U(:).^(model.flux_pdeg-1);   
  flux_lin.Vy = RV(2,:)' .* model.flux_pdeg .* U(:).^(model.flux_pdeg-1);
  %  flux.lambda = 1/model.vmax;
  
else % perform forward difference
  epsilon = max(abs(U)) * 1e-2;
  if epsilon < eps; % U completely 0
    epsilon = 1e-2;
  end;
  
  fluxU = conv_flux(model,U,X,Y);
  fluxU2 = conv_flux(model,U+epsilon,X,Y);
  
  flux_lin.Vx = (fluxU2.Fx-fluxU.Fx)/epsilon;
  flux_lin.Vy = (fluxU2.Fy-fluxU.Fy)/epsilon;
  flux_lin.lambda = fluxU.lambda;
end;

if decomp_mode>0 & respected_decomp_mode==0
  error('function does not support affine decomposition!');
end;
 

%| \docupdate 
