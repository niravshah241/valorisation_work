function Udirichlet = dirichlet_values(model,X,Y)
%UDIRICHLET = DIRICHLET_VALUES([X],[Y], MODEL)
% Examples
%  dirichlet_values([0,1,2],[1,1,1],struct('name_dirichlet_values', 'homogeneous', ...
%                                        'c_dir', 1))
%  dirichlet_values([0:0.1:1],[0],struct('name_dirichlet_values', 'xstripes', ...
%                                        'c_dir', [0 1 2], ...
%                                        'dir_borders', [0.3 0.6]))
%  dirichlet_values([0:0.1:1],[0],struct('name_dirichlet_values', 'box', ...
%                                        'c_dir', 1, ...
%                                        'dir_box_xrange', [0.3 0.6], ...
%                                        'dir_box_yrange', [-0.1 0.1]))
%
% function computing dirichlet-values by pointwise evaluations  
% 'zero' Udirichlet = 0
% 'homogeneous' Udirichlet = c_dir
% 'leftright' Udirichlet = c_dir_left for x <= dir_middle 
%             Udirichlet = c_dir_right for x > dir_middle   
% 'xstripes' values constant in stripes divided by x-coordinats
%            n+1 values are given in c_dir
%            the separating coodrinates in dir_borders
% 'box': values constant model.c_dir in box given by
%         model.dir_box_xrange and dir_box_yrange
% 'weighted_boxes': values constant model.c_dir in box given by
%       (beta=1):  model.dir_box_xrange{1} and dir_box_yrange{1} 
%       (beta=0):  model.dir_box_xrange{2} and dir_box_yrange{2} 
%            for intermediate betas, the two boxes are weighted
% 'gauss_convcomb': convex combination of gaussian distributions at 
%                   left border
% 'function_ptr' : dirichlet function given by function
%                  pointer to complete evaluation:
%                      model.dirichlet_values_ptr
%                  Function must
%                  allow call U = f(X,Y,model), where X, Y are
%                  vectors of the same size. Return is a vector of
%                  values in the corresponding points.
% 'decomp_function_ptr' : dirichlet function given by function
%                  pointer to components and coefficients: 
%                      model.dirichlet_values_coefficients_ptr 
%                  Function must
%                  allow call U = f([],[],model)
%                  Return is a vector of length Q with (parameter
%                  dependent) coefficients in U.
%                      model.dirichlet_values_components_ptr 
%                  Functions must
%                  allow call U = f(X,Y,model), where X, Y are
%                  vectors of the same size. 
%                  Return of components is a cell array of vectors of
%                  the same size as X and Y with the point values
%                  of the components. Linear combination of
%                  components by coefficients then yields the
%                  complete evaluation.
%  
% required fields of model:
%   name_dirichlet_values : 'zero', 'homogeneous', 'leftright', 'uplow'
%   c_dir(s)              : boundary value(s) in case of homogeneous dirichlet
%                           value or stripes or box
%   c_dir_left            : left dirichlet value
%   c_dir_right           : left dirichlet value
%   c_dir_up              : upper dirichlet value
%   c_dir_low             : lower dirichlet value
%   dir_middle            : coordinate for separatin between two domains
%   dir_borders           : coordinates for separation between stripes
%   dir_box_xrange        : x-coordinate interval for dirichlet box
%   dir_box_yrange        : y-coordinate interval for dirichlet box
%
% in case of gld-dirichlet-value: upper wall value c_dir, remaining 0.0  
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
    
  if isequal(model.name_dirichlet_values,'zero')
    if decomp_mode == 2
      Udirichlet = 0;
    elseif decomp_mode == 1      
      Udirichlet = {[]};
    elseif decomp_mode == 0      
      Udirichlet = zeros(length(X),1);
    else
      error('unknown decomp_mode');
    end;
    respected_decomp_mode = 1;
  elseif isequal(model.name_dirichlet_values,'homogeneous')
    Udirichlet = model.c_dir * ones(length(X),1);
  elseif isequal(model.name_dirichlet_values,'leftright')
    i = (X <= model.dir_middle);
    U0 = i;
    U1 = (1-i);
    if decomp_mode == 2
      Udirichlet = [model.c_dir_left, model.c_dir_right];
    elseif decomp_mode == 1      
      Udirichlet = {U0,U1};
    elseif decomp_mode == 0      
      Udirichlet = U0 * model.c_dir_left  + U1 * model.c_dir_right;
    end;    
    respected_decomp_mode = 1;
  elseif isequal(model.name_dirichlet_values,'uplow')
    i = (Y <= model.dir_middle);
    if decomp_mode == 2
      Udirichlet = [model.c_dir_low, model.c_dir_up];
    elseif decomp_mode == 1
      Udirichlet = {i, (1-i)};
    elseif decomp_mode == 0
      Udirichlet = model.c_dir_low * i + (1-i) * model.c_dir_up;
    end;
    respected_decomp_mode = 1;
  elseif isequal(model.name_dirichlet_values,'xstripes')
    % setzen aller Streifen von links nach rechts, (mehrmals ueberschrieben) 
    Udirichlet = ones(size(X)) * model.c_dir(1); 
    for n = 2:length(model.c_dir)
      fi = X > model.dir_borders(n-1);
      Udirichlet(fi) = model.c_dir(n);
    end;
  elseif isequal(model.name_dirichlet_values,'box')
    Udirichlet = zeros(size(X)); 
    fi = X > model.dir_box_xrange(1) & ...
	       X < model.dir_box_xrange(2) & ...
	       Y > model.dir_box_yrange(1) & ...
	       Y < model.dir_box_yrange(2);
    Udirichlet(fi) = model.c_dir;
  elseif isequal(model.name_dirichlet_values,'weighted_boxes')
    if decomp_mode == 0
      Udirichlet = zeros(size(X)); 
      xrange = model.dir_box_xrange{1};
      yrange = model.dir_box_yrange{1};
      fi = X > xrange(1) & ...
        X < xrange(2) & ...
        Y > yrange(1) & ...
        Y < yrange(2);
      Udirichlet(fi) = model.c_dir*model.beta;
      xrange = model.dir_box_xrange{2};
      yrange = model.dir_box_yrange{2};
      fi = X > xrange(1) & ...
        X < xrange(2) & ...
        Y > yrange(1) & ...
        Y < yrange(2);
      Udirichlet(fi) = model.c_dir*(1-model.beta);
    elseif decomp_mode == 1
      % two components if beta is in mu, otherwise one component
      Udirichlet = cell(2,1);
      for q=1:2
        Udirichlet{q} = zeros(size(X)); 
        xrange = model.dir_box_xrange{q};
        yrange = model.dir_box_yrange{q};
        fi = X > xrange(1) & ...
          X < xrange(2) & ...
          Y > yrange(1) & ...
          Y < yrange(2);
        Udirichlet{q}(fi) = 1;
      end;
      if ~ismember('beta',model.mu_names) 
        % merge to one component
        Udirichlet = {model.beta*Udirichlet{1} + ...
		      (1-model.beta)*Udirichlet{2}};
      end;
    else % decomp_mode = 2
      if ~ismember('beta',model.mu_names) 
        Udirichlet = model.c_dir;
      else
        Udirichlet = model.c_dir * [model.beta; 1-model.beta];
      end;
    end;
    respected_decomp_mode = 1;
    
  elseif isequal(model.name_dirichlet_values,'gauss_convcomb')    
    xscale = model.udir_xscale;
    if decomp_mode == 0
      a = model.udir_amplitude;
      gauss = cell(6);
      for i = 1:6
	gauss{i} = a * exp(-100 * ((Y-0.1*i).^2 + (X*xscale).^2));
      end;
      
      mu4 = model.udir_height;
      if mu4>=1 &&  mu4<=2
	coeff = [2-mu4, mu4-1, 0, 0, 0, 0];
      elseif mu4>2 &&  mu4<=3
	coeff = [0, 3-mu4, mu4-2, 0, 0, 0];
      elseif mu4>3 &&  mu4<=4
	coeff = [0, 0, 4-mu4, mu4-3, 0, 0];
      elseif mu4>4 &&  mu4<=5
	coeff = [0, 0, 0, 5-mu4, mu4-4, 0];
      elseif mu4>5 &&  mu4<=6	
	coeff = [0, 0, 0, 0, 6-mu4, mu4-5];
      end;
      
      Udirichlet  = zeros(size(X));
      for i = 1:6
	Udirichlet = Udirichlet + gauss{i} * coeff(i);
      end;
    elseif decomp_mode == 1
      a = model.udir_amplitude;
      Udirichlet = cell(6);
      for i = 1:6
	Udirichlet{i} = a * exp(-100 * ((Y-0.1*i).^2 + (X*xscale).^2));
      end;
    elseif decomp_mode == 2 % coefficients
      mu4 = model.udir_height;
      if mu4>=1 &&  mu4<=2
	Udirichlet = [2-mu4, mu4-1, 0, 0, 0, 0];
      elseif mu4>2 &&  mu4<=3
	Udirichlet = [0, 3-mu4, mu4-2, 0, 0, 0];
      elseif mu4>3 &&  mu4<=4
	Udirichlet = [0, 0, 4-mu4, mu4-3, 0, 0];
      elseif mu4>4 &&  mu4<=5
	Udirichlet = [0, 0, 0, 5-mu4, mu4-4, 0];
      elseif mu4>5 &&  mu4<=6	
	Udirichlet = [0, 0, 0, 0, 6-mu4, mu4-5];
      else
	error('udir_height out of range!');
      end;
    else
      error('dirichlet value not implemented for complete or components!');
    end;
    respected_decomp_mode = 1;

  elseif isequal(model.name_dirichlet_values,'decomp_function_ptr')
    if decomp_mode == 2
      Udirichlet = model.dirichlet_values_coefficients_ptr([],[],model);
    elseif decomp_mode == 1; % components
      Udirichlet = model.dirichlet_values_components_ptr(X,Y,model);
    else % decomp_mode = 0, complete
      Ucoefficients = model.dirichlet_values_coefficients_ptr([],[],model);
      Ucomponents = model.dirichlet_values_components_ptr(X,Y,model); 
      Udirichlet = lincomb_sequence(Ucomponents,Ucoefficients);	  
    end;    
    respected_decomp_mode = 1;
    
  elseif isequal(model.name_dirichlet_values,'function_ptr')
    Udirichlet = model.dirichlet_values_ptr(X,Y,model);
    respected_decomp_mode = 0;
    
  else
    error('unknown name_dirichlet_values');
  end;

  if decomp_mode>0 && respected_decomp_mode==0
    error('function does not support affine decomposition!');
  end;
 
%| \docupdate 
