function U0 = init_values_old_to_be_split_in_files(X,Y, params)
%function U0 = init_values([X],[Y], params)
%
% function constructing the initial values of the convection diffusion
% problem in the specified points and parameters
%
% 'blobs':   u_0 = beta*exp(-gamma*((x-0.25).^2+(y-0.5).^2))* Xi(B_r) ...
%                  + (1-beta)*exp(-gamma*((x-0.25).^2+(y-0.7).^2))* Xi(B_r)
% 'homogeneous' u_0 = beta
% 'as_dirichlet' : evaluate dirichlet-data-function at time 0  
% 'waveproduct'  : product of two sinus waves in both coordinate directions
% 'function_ptr' : init function given by function
%                  pointer to complete evaluation:
%                      params.init_values_ptr
%                  Function must
%                  allow call U = f(X,Y,params), where X, Y are
%                  vectors of the same size. Return is a vector of
%                  values in the corresponding points.
% 'decomp_function_ptr' : init function given by function
%                  pointer to components and coefficients: 
%                      params.init_values_coefficients_ptr 
%                  Function must
%                  allow call U = f([],[],params)
%                  Return is a vector of length Q with (parameter
%                  dependent) coefficients in U.
%                      params.init_values_components_ptr 
%                  Functions must
%                  allow call U = f(X,Y,params), where X, Y are
%                  vectors of the same size. 
%                  Return of components is a cell array of vectors of
%                  the same size as X and Y with the point values
%                  of the components. Linear combination of
%                  components by coefficients then yields the
%                  complete evaluation.
%
% required fields of params:
% name_init_values     : 'blobs', 'homogeneous', 'wave', 'leftright', 
%                        'as_dirichlet'
%
% if name_init_values == 'gradient_box'
%   c_init_up       : constant value on upper border of box
%                     (unit square)
%   c_init_lo       : constant value on lower border of box
%
% if name_init_values == 'blobs'
% radius             : radius of blob-cutoff circle
% gamma              : spread of initial-value gaussian
% beta               : value between 0 and 1 weighting two gauss-blobs
%
% if name_init_values == 'homogeneous'
% c_init             : constant value in field
%
% if name_init_values == 'wave'
% a wave of values between 0 and c_init 
%           c_init* 0.5 * (sin(freq_x * x + freq_y * y) + 1)
% c_init             : maximum value in field
% freq_x            : frequency in x-direction
% freq_y            : frequency in y-direction
%
% if name_init_values == 'waveproduct'
% a product of axis-dependent waves of values between 
%           c_init_min and c_init_max 
% c_init_max             : maximum value in field
% c_init_min             : minimum value in field
% c_init_freq_x           : frequency in x-direction
% c_init_freq_y           : frequency in y-direction
% c_init_phase_x          : phase shift
% c_init_phase_y          : phase shift
%
% if name_init_values == 'leftright'
% c_init_left      : constant value in field left of border
% c_init_right     : constant value in field left of border
% c_init_middle    : x-coordinate separation  
%
% if name_init_values == 'as_dirichlet'
%          parameters as required by dirichlet function
%
% Function supports affine decomposition, i.e. different operation modes
% guided by optional field decomp_mode in params. See also the
% contents.txt for general explanation
%
% optional fields of params:
%   mu_names : names of fields to be regarded as parameters in vector mu
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


% Bernard Haasdonk 5.10.2005

% determine affine_decomposition_mode as integer
decomp_mode = params.decomp_mode;
% flag indicating whether the computation respected the decomposition
respected_decomp_mode = 0;

if isequal(params.name_init_values,'blobs')
  B1 = ((X(:)-0.25).^2+(Y(:)-0.5).^2);
  B2 = ((X(:)-0.25).^2+(Y(:)-0.7).^2);
  I1 = (B1<params.radius.^2);
  I2 = (B2<params.radius.^2);
  if decomp_mode == 0 
    U0 = params.beta*exp(-params.gamma*B1).*I1 ...
	 + (1-params.beta)* exp(-params.gamma*B2).*I2;
  elseif decomp_mode == 1 
    if ismember('gamma',params.mu_names) | ...
	  ismember('radius',params.mu_names)
      error('affine decomp with respect to mu_names not possible!');
    end;
    % two components
    U0{1} = exp(-params.gamma*B1).*I1;
    U0{2} = exp(-params.gamma*B2).*I2;      
  else % decomp_mode == 2 
    U0 = [params.beta, 1-params.beta];
  end;  
  respected_decomp_mode = 1;
  
elseif isequal(params.name_init_values,'leftright')
  i = (X <= params.init_middle);
  U0 = params.c_init_left * i + (1-i) * params.c_init_right;
  
elseif isequal(params.name_init_values,'decomp_function_ptr')
  if decomp_mode == 2
    U0 = params.init_values_coefficients_ptr([],[],params);
  elseif decomp_mode == 1; % components
    U0 = params.init_values_components_ptr(X,Y,params);
  else % decomp_mode = 0, complete
    Ucoefficients = params.init_values_coefficients_ptr([],[],params);
    Ucomponents = params.init_values_components_ptr(X,Y,params); 
    U0 = lincomb_sequence(Ucomponents,Ucoefficients);	  
  end;    
  respected_decomp_mode = 1;
  
elseif isequal(params.name_init_values,'function_ptr')
  U0 = params.init_values_ptr(X,Y,params);
  respected_decomp_mode = 0;  

else
  error('unknown name_init_values');
end;

if decomp_mode>0 && respected_decomp_mode==0
  error('function does not support affine decomposition!');
end;

% TO BE ADJUSTED TO NEW SYNTAX
%| \docupdate 
