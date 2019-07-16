function U0 = init_values(X,Y, params)
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
% 'grey_image'   : load grey value image as initial data
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
% if name_init_values == 'grey_image'
%          a grey value image is loaded mapping 0 to 0 and 255 to c_init
%          with bilinear interpolation 
% c_init             : parameter weighting the components
% c_init_filename    : name of image file
% c_init_xmin        : image is mapped to specified rectangle 
% c_init_xmax        : image is mapped to specified rectangle 
% c_init_ymin        : image is mapped to specified rectangle 
% c_init_ymax        : image is mapped to specified rectangle 
% c_init_xdivisions  : vector of divisions, i.e intervals,
%                      e.g. borders of letters in an image
%                      determines number of components. c_init = 0
%                      => full weight to leftmost part of image
%                      c_init = 1 => full weight to complete image
%
% Function supports affine decomposition, i.e. different operation modes
% guided by optional field affine_decomp_mode in params. See also the 
% contents.txt for general explanation
%
% optional fields of params:
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

decomp_mode = params.decomp_mode;
% flag indicating whether the computation respected the decomposition
respected_decomp_mode = 0;

if isequal(params.name_init_values,'gradient_box')
  % definition of the box
  width = 0.2;
  xmid = 0.45;
  ymid = 0.55;
  height = 0.7;
  ymin = ymid-height/2;
  ymax = ymid+ height/2;
  if decomp_mode == 0 
    BX = abs(X(:)-xmid);
    BY = abs(Y(:)-ymid);
    I = (BX<width/2) & (BY<height/2);
    % y = ymin =>  0 + c_init_lo
    % y = ymax =>  c_init_up + 0
    U0 = (params.c_init_up * (Y(:)-ymin) + ...
	  + params.c_init_lo * (ymax-Y(:))) .* I / (ymax-ymin);
  elseif decomp_mode == 1 
    % two components
    BX = abs(X(:)-xmid);
    BY = abs(Y(:)-ymid);
    I = (BX<width/2) & (BY<height/2);
    U0{1} = (Y(:)-ymin) .* I / (ymax-ymin);
    U0{2} = (ymax-Y(:)) .* I / (ymax-ymin);
  else % decomp_mode == 2 
    U0 = [params.c_init_up, params.c_init_lo];
  end;  
  respected_decomp_mode = 1;
elseif isequal(params.name_init_values,'blobs')
  B1 = ((X(:)-0.25).^2+(Y(:)-0.5).^2);
  B2 = ((X(:)-0.25).^2+(Y(:)-0.7).^2);
  I1 = (B1<params.radius.^2);
  I2 = (B2<params.radius.^2);
  if decomp_mode == 0 
    U0 = params.beta*exp(-params.gamma*B1).*I1 ...
	 + (1-params.beta)* exp(-params.gamma*B2).*I2;
  elseif decomp_mode == 1 
    if ismember('gamma',params.mu_names) || ...
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
elseif isequal(params.name_init_values,'homogeneous')
    if decomp_mode == 0
        U0 = params.c_init * ones(length(X),1);
    elseif decomp_mode == 1
        U0 = { ones(length(X),1) };
    else % decomp_mode == 2
        U0 = params.c_init;
    end;
    respected_decomp_mode = 1;
elseif isequal(params.name_init_values,'wave')
  if decomp_mode == 0
    U0 = params.c_init * 0.5 * (sin(params.freq_x * X(:) + ...
				    params.freq_y * Y(:)  ) + 1 );
  elseif decomp_mode == 1 
    % single component independent of mu_names
    if ismember('freq_x',params.mu_names) || ...
	  ismember('freq.y',params.mu_names)
      error('affine decomp with respect to mu_names not possible!');
    end;
    U0 = {0.5 * (sin(params.freq_x * X(:) + ...
		     params.freq_y * Y(:)  ) + 1 )};      
  else % decomp_mode == 2 
    U0 = params.c_init;
  end;
  respected_decomp_mode = 1;
  
  elseif isequal(params.name_init_values,'transformed_blobs')
%    [X_trans, Y_trans] = geometry_transformation(X,Y,params);
%    B = cell(1,5);
%    I = cell(1,5);
%    B{1} = ((X_trans(:) - 0.5).^2 +(Y_trans(:)-0.7).^2);
    if any(ismember(fieldnames(params), 'blob_height'))
      h = params.blob_height;
    else
      h = 0.15;
    end
    B{1} = ((X(:) - 0.5).^2 +(Y(:)-0.7).^2);
    B{2} = ((X(:) - 0.8).^2 +(Y(:)-0.9).^2);
    B{3} = ((X(:) - 0.8).^2 +(Y(:)-0.2).^2);
    B{4} = ((X(:) - 0.2).^2 +(Y(:)-0.2).^2);
    B{5} = ((X(:) - 0.2).^2 +(Y(:)-0.9).^2);
    B{6} = ((X(:) - 0.4).^2 +(Y(:)-0.4).^2);
    B{7} = ((X(:) - 0.6).^2 +(Y(:)-0.6).^2);
                    
    I=cell(1,7);
    % p(mu) = mu*(1-x)
    % p_mu = spline_select(params);
    ztrans = Y.*(1-X);
    for i=1:7
        I{i} = (B{i} < 0.05.^2);
    end
    if decomp_mode == 0
      U0 = h .* I{1} + params.c_init * ones(length(X),1) + params.gravity * (Y + params.hill_height*ztrans);
      for i = 2:7
        U0 = U0 + h .* exp(-B{i}) .* I{i};
      end
    elseif decomp_mode == 1
      U0 = cell(1,8);
      U0{1} = I{1};
        for i = 2:7
        U0{i} = I{i};
        end
      U0{8} = ones(length(X),1);
      if params.gravity ~= 0
        U0{9} = Y;
        U0{10} = ztrans;
      end
    else
      U0 = [ h, h, h, h, h, h, h, params.c_init];
      U0 = [ U0, params.gravity, params.gravity * params.hill_height ];
    end
    respected_decomp_mode = 1;
elseif isequal(params.name_init_values,'waveproduct')
  if (decomp_mode == 0) || (decomp_mode ==1)
    Uinit = params.c_init_min + ...
	    (params.c_init_max-params.c_init_min) * ...
	    0.5 * (sin_sym(params.c_init_freq_x * X(:) +  params.c_init_phase_x) ...
		   .* sin_sym(params.c_init_freq_y * Y(:) +  params.c_init_phase_y) ...
		   + 1);
  end;
%  if ismember('c_init_freq_x',params.mu_names) | ...
%	ismember('c_init_freq_y',params.mu_names) | ...
%	ismember('c_init_min',params.mu_names) | ...
%	ismember('c_init_max',params.mu_names) | ...
%	ismember('c_init_phase_x',params.mu_names) | ...
%	ismember('c_init_phase_y',params.mu_names) 
%    error('affine decomp with respect to mu_names not possible!');
%  end;
  switch decomp_mode 
   case 0
    U0 = Uinit;
   case 1
    U0 = {Uinit};
   case 2
    U0 = 1;
  end;
  respected_decomp_mode = 1;
  
elseif isequal(params.name_init_values,'leftright')
  i = (X <= params.init_middle);
  U0 = params.c_init_left * i + (1-i) * params.c_init_right;
elseif isequal(params.name_init_values,'as_dirichlet')
  U0 = dirichlet_values(X,Y,params);
  % if as_dirichlet not decomposable => error there and not here.
  respected_decomp_mode = 1;

elseif isequal(params.name_init_values,'grey_image')
  ncomp = 1+length(params.c_init_xdivisions);
  if decomp_mode ~= 2
    % determine components and perform linear combination
    img = imread(params.c_init_filename)/255;
    %erstmal ohne bilineare Interpolation:
    XI = (X(:)-params.c_init_xmin)/(params.c_init_xmax- ...
				    params.c_init_xmin)* size(img,2);
    XI = max(XI,0);
    XI = min(XI,size(img,2));
    XI = ceil(XI);
    fi = find(XI==0);
    if ~isempty(fi)
      XI(fi)= 1;
    end;    
    YI = (Y(:)-params.c_init_ymin)/(params.c_init_ymax- ...
				    params.c_init_ymin)* size(img,1);
    YI = max(YI,0);
    YI = min(YI,size(img,1));
    % reflect y coordinate
    YI = size(img,1)-YI;
    YI = ceil(YI);
    fi = find(YI==0);
    if ~isempty(fi)
      YI(fi)= 1;
    end;
    
    ind = sub2ind(size(img),YI,XI);
    I = img(ind);
    
    U0comp = cell(ncomp,1);
    xdivs = [0;params.c_init_xdivisions(:);params.c_init_xmax];
    for i = 1:ncomp
      U0 = zeros(length(X(:)),1);
      fi = find( (X(:)<xdivs(i+1)) & (X(:)>=xdivs(i)));
      if ~isempty(fi)
	U0(fi) = I(fi);
      end;
      U0comp{i} = U0;
    end;
  end;
    
  % determine coefficients:
  %coeff = zeros(ncomp,1);
  coeff = (0:(ncomp-1))/ncomp;
  coeff = (params.c_init - coeff)*ncomp;
  coeff = min(coeff,1);
  coeff = max(coeff,0);
  
  % now return the correct quantity
  if decomp_mode == 2 
    U0 = coeff;
  elseif decomp_mode == 0 
    U0 = lincomb_sequence(U0comp,coeff);
  elseif decomp_mode == 1 
    U0 = U0comp;
  end;
  respected_decomp_mode = 1;
  
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
