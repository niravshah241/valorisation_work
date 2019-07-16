function FNneu = neuman_values(model,Xneu,Yneu,Uneu,Nxneu,Nyneu)
%function FNneu = neuman_values(model,[Xneu],[Yneu],[Uneu],[Nxneu],[Nyneu])
%
% function computing neuman-values by pointwise evaluation. 
%
% Xneu, Yneu: coordinate vectors of points to be evaluated
% Uneu:       corresponding U-value in case of a u-dependent flux
% Nxneu,Nyneu: corresponding unit normal vectors
% 
% required field of model:  name_neuman_values:   
% 'zero'  FNneu = 0, arguments except Xneu may be ignored
% 'homogeneous' FNneu = c_neu, arguments except Xneu may be ignored
% 'outflow': FNneu = F(Uneu, Xneu, Yneu) * n
% 'pressure_gdl' FNeu = no_flow on upper and lower boundary. left and
%                right is a parabolic velocity profile with maximum c_neu_max.
%                Uneu and Nxneu and Nyneu may be ignored.
% 'rightflow': outflow to right, noflow to upper and lower
%
% required fields of model:
%   name_neuman_values  : 'zero', 'homogeneous', 'outflow', 'pressure_gdl'
%   c_neu :               boundary value in case of homogeneous 
%                         dirichlet value
%   c_neu_max : maximum flow value for 'pressure_gdl'
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


% Bernard Haasdonk 11.4.2006

% determine affine_decomposition_mode as integer  
decomp_mode = model.decomp_mode;
%decomp_mode = get_affine_decomp_mode(model);
% flag indicating whether the computation respected the decomposition
respected_decomp_mode = 0;

  if isequal(model.name_neuman_values,'zero')
    if decomp_mode == 2
      % cannot detect, whether any points will be available later,
      % so at least one component must be available.
      FNneu = 0;
    elseif decomp_mode == 1
      % cannot detect in mode 2, whether any points will be available later,
      % so at least one component must be available.      
      FNneu = {[]};
    elseif decomp_mode == 0
      FNneu = zeros(size(Xneu(:)));
    else
      error('unknown decomposition mode');
    end;
    respected_decomp_mode = 1;
  elseif isequal(model.name_neuman_values,'homogeneous')
    FNneu = model.c_neu * ones(size(Xneu(:)));
  elseif isequal(model.name_neuman_values,'outflow')
    Fneu = conv_flux(model,Uneu, Xneu, Yneu);
    if decomp_mode == 0
      FNneu = Fneu.Fx .* Nxneu(:) + Fneu.Fy .* Nyneu(:);
    elseif decomp_mode == 1
      Q_Fneu = length(Fneu);
      FNneu = cell(Q_Fneu,1);
      FNneu{:} = zeros(size(Xneu(:)));
      for q = 1:Q_Fneu
	FNneu{q} =  Fneu{q}.Fx.*Nxneu(:) + Fneu{q}.Fy.*Nyneu(:);
      end;
    else % mode 2
      % error('Not sure what to do in this case.')
      FNneu = Fneu;  % perhaps this is right?
    end;
    respected_decomp_mode = 1;    
  elseif isequal(model.name_neuman_values,'pressure_gdl')
    % set all to zero
    FNneu = zeros(size(Xneu(:))); 			
    % left boundary negative velocity -> positive pressure gradient      
    i = find(Xneu<model.xrange(1)+eps);
    FNneu(i) = model.c_neu_max * ...
	(model.yrange(2)-model.yrange(1))^(-2) * ...
	(Yneu(i)-model.yrange(1)).*(model.yrange(2)-Yneu(i))*4; 
    % right boundary positive velocity -> negative pressure gradient    
    i = find(Xneu>model.xrange(2)-eps);
    FNneu(i) = -model.c_neu_max * ...
	(model.yrange(2)-model.yrange(1))^(-2) * ...
	(Yneu(i)-model.yrange(1)).*(model.yrange(2)-Yneu(i))*4; 
  elseif isequal(model.name_neuman_values,'rightflow')
    if decomp_mode < 2
      i = find(Xneu>model.xrange(2)-eps);   

      % in case of filecaching == 2, the correct velocity file must
      % be generated!!
      if isfield(model,'filecache_velocity_matrixfile_extract') & ...
	    (model.filecache_velocity_matrixfile_extract == 2);
	model.velocity_matrixfile = ... 
	    cache_velocity_matrixfile_extract(...
		Xneu(i),Yneu(i),'outflow', model);
      end;
      Fneu = conv_flux(model,Uneu(i), Xneu(i), Yneu(i));

    end;
    if decomp_mode == 0 % none
                  
      FNneu = zeros(size(Xneu(:))); 		
      FNneu(i) = Fneu.Fx.*Nxneu(i) + Fneu.Fy.*Nyneu(i);
    elseif decomp_mode == 1 % components
      Q_Fneu = length(Fneu);
      FNneu = cell(Q_Fneu,1);
      FNneu{:} = zeros(size(Xneu(:)));
      for q = 1:Q_Fneu
	FNneu{q}(i) =  Fneu{q}.Fx.*Nxneu(i) + Fneu{q}.Fy.*Nyneu(i);
      end;
      % check dependency on xrange
      if ismember('xrange',model.mu_names)
	error('affine decomp with respect to mu_names not possible!');
      end;      
    else % decomp_mode == 2 -> coefficients
      Fneu = conv_flux(model,[], [], []);
      FNneu = Fneu; % simple identical coefficients!
    end;
    respected_decomp_mode = 1;
  else
    error('unknown name_neuman_values');
  end;
  
  if decomp_mode>0 & respected_decomp_mode==0
    error('function does not support affine decomposition!');
  end;
  
%| \docupdate 
