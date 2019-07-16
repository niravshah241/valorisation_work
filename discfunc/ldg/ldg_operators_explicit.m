function [L_I,L_E,b] = ldg_operators_explicit(model,model_data)
%function [L_I,L_E,b] = ldg_operators_explicit(model,model_data)
% 
% function computing the time evolution matrices for a ldg time
% step  L_I * Unew = L_E * U + b
% result are two sparse matrices L_I, L_E and an offset vector b
% Function supports affine decomposition, i.e. different operation modes
% guided by optional field decomp_mode in model. See also the 
% contents.txt for general explanation
%
% required fields of model:
%   dt      : timestep
%   t       : current time, i.e. discretization is taking place in t and
%             t+dt for explicit and implicit parts respectively
%   verbose : integer indicating the verbosity level
%
% optional fields of model:
%   mu_names : names of fields to be regarded as parameters in vector mu

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


% Bernard Haasdonk 27.8.2009
  
% determine affine_decomposition_mode as integer  

disp('to be adjusted!');
keyboard;

decomp_mode = model.decomp_mode;

% call operator function with or without grid
if decomp_mode < 2 
  n = length(model_data.grid.A);
end;
%  [L_I_diff, bdir_I_diff] = fv_operators_diff_implicit(model,model_data);
%  [L_E_diff, bdir_E_diff] = fv_operators_diff_explicit(model,model_data);
[L_E_conv, bdir_E_conv] = ldg_operators_adv_explicit(model, ...
						  model_data);
%  [L_I_conv, bdir_I_conv] = ldg_operators_conv_implicit(model,model_data);
%  [L_E_neu,  bneu_E]      = fv_operators_neuman_explicit(model,model_data);
%  [L_I_neu,  bneu_I]      = fv_operators_neuman_implicit(model,model_data);

% check that offline and online simulation generate identical
% number of components, for this run routine twice and compare numbers:
%if model.verbose>=10 && decomp_mode >=1
%  vnames = {'L_I_diff','bdir_I_diff','L_E_diff','bdir_E_diff',...
%	   'L_E_conv','bdir_E_conv','L_I_conv','bdir_I_conv',...
%	   'L_E_neu','bneu_E','L_E_neu','bneu_I'};
%  for v = 1:length(vnames)
%    disp([vnames{v},': ',num2str(size(eval(vnames{v})))]);
%  end
%  keyboard;
%end;

%if decomp_mode==0  
%  % check for condition: 1-dt * L_E_conv(i,i) >= 0 as this is the 
%  % coefficient of u_i^k in the convex-combination representation of u_i^(k+1) 
%  % maximum possible dt_max = min ( 1/L_E_conv(i,i))
%  dt_max = min(full(diag(L_E_conv)).^(-1));
%  if model.verbose>=9
%    if model.dt > dt_max
%      disp(['current dt = ',num2str(model.dt),...
%	    ' is larger than cfl-condition!']);
%      disp(['cfl: dt_max = ',num2str(dt_max)]);
%    end;
%    if model.dt < dt_max * 0.1;
%      disp(['current dt = ',num2str(model.dt),' can be chosen much larger!']);
%      disp(['cfl: dt_max = ',num2str(dt_max)]);
%    end;
%  end;
%end;

if decomp_mode == 2
  % decomp_mode == 2 -> sigmas are required
  %  L_I = [1; model.dt* L_I_diff(:); ...
  %	 model.dt * L_I_conv(:); model.dt * L_I_neu(:)];
  %  L_E = [1; - model.dt* L_E_diff(:); ...
  %	 -model.dt * L_E_conv(:); -model.dt * L_E_neu(:)];
  %  b = model.dt * [ bdir_E_diff(:) ; bdir_I_diff(:); ...
  %		    bdir_E_conv(:) ; bdir_I_conv(:); ...
  %		    bneu_E(:); bneu_I(:)];
  L_I = 1;
  L_E = [1; -model.dt * L_E_conv(:)];
  b = model.dt * bdir_E_conv(:);
elseif decomp_mode == 1
  % in case of components: concatenate eye and the other lists
  %  L_I = [{speye(n)}; L_I_diff(:); L_I_conv(:); L_I_neu(:)];
  %  L_E = [{speye(n)}; L_E_diff(:); L_E_conv(:); L_E_neu(:)];
  %  b = [bdir_E_diff(:); bdir_I_diff(:); ...
  %       bdir_E_conv(:); bdir_I_conv(:); ...
  %       bneu_E(:)     ; bneu_I(:)];
  L_I = {speye(n)}; 
  L_E = [{speye(n)}; L_E_conv(:)];
  b = bdir_E_conv(:);
  %  keyboard;
else % decomp_mode = 0 => complete
  L_I = speye(n);
  L_E = speye(n) - model.dt * L_E_conv;
  b  = model.dt * (bdir_E_conv);
end;

%if decomp_mode == 0
%  % check L_E contributions:
%  % matrix L should have positive diagonal and nonpositive off-diagonal
%  % with row-sum between 0 and 1
%  L = L_E_diff + L_E_conv + L_E_neu;
%  
%  if model.verbose>=10
%    if ~isempty(find(diag(L)<0))
%      disp('warning: explicit matrix has negative diagonal entries!!');
%      if model.verbose>=10
%	keyboard;
%      end;
%    end;
%    Lo = L-diag(diag(L));
%    if ~isempty(find(Lo>0))
%      disp('warning: explicit matrix has positive off-diagona entries!!');
%      if model.verbose>=10
%	keyboard;
%      end;
%    end;
%    if ~isempty(find(sum(L,2)<0))
%      disp('warning: explicit matrix has negative row sum contributions!!');
%      if model.verbose>=10
%	keyboard;
%      end;
%    end;
%  end;
%  
% % only report this time consuming computation once:
%  if (model.t==0) && (model.verbose >= 10)
%    fv_estimate_operator_norms(model,model_data);
%  end;
%end;

%| \docupdate 
