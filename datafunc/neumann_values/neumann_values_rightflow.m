function FNneu = neumann_values_rightflow(glob,U,normals,params)
%function FNneu = neumann_values_rightflow(glob,U,normals,params)
%
% function computing neuman-values by pointwise evaluation. 
%
% glob: columnwise coordinate vectors of global points to be evaluated
% Uneu: columnwise U values in points in case of a u-dependent flux
% normals: Columnwise corresponding unit normal vectors
% 
% 'rightflow': outflow to right, noflow to upper and lower
%
% Function supports affine decomposition, i.e. different operation modes
% guided by optional field decomp_mode in params. See also the
% contents.txt for general explanation

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


% Bernard Haasdonk 3.9.2009

% determine affine_decomposition_mode as integer  
% glob column check
% normals column check
if params.debug
  if ~isempty(normals) && size(glob,1) < size(glob,2)
    warning('coordinates in variable normals are given row-wise, but expected them to be column-wise');
    if params.debug > 2
      keyboard;
    end
  end
end
if params.debug
  if ~isempty(glob) && size(glob,1) < size(glob,2)
    warning('coordinates in variable glob are given row-wise, but expected them to be column-wise');
    if params.debug > 2
      keyboard;
    end
  end
end
decomp_mode = params.decomp_mode;

if decomp_mode < 2
  i = find(glob(:,1)>params.xrange(2)-eps);   
  
  % in case of filecaching == 2, the correct velocity file must
  % be generated!!
  if params.filecache_velocity_matrixfile_extract == 2;
    params.velocity_matrixfile = ... 
	cache_velocity_matrixfile_extract(...
	    glob(i,1),glob(i,2),'outflow', params);
  end;
  Fneu = params.conv_flux_ptr(glob(i,:),U(i,:),params);
  
end;
if decomp_mode == 0 % none
  
  FNneu = zeros(size(glob,1),1); 		
  FNneu(i) = Fneu(:,1).*normals(i,1) + Fneu(:,2).*normals(i,2);
elseif decomp_mode == 1 % components
  Q_Fneu = length(Fneu);
  FNneu = cell(Q_Fneu,1);
  FNneu{:} = zeros(size(glob,1),1);
  for q = 1:Q_Fneu
    FNneu{q}(i) =  Fneu{q}(:,1).*normals(i,1) + Fneu{q}(:,2).*normals(i,2);
  end;
  % check dependency on xrange
  if ismember('xrange',params.mu_names)
    error('affine decomp with respect to mu_names not possible!');
  end;      
else % decomp_mode == 2 -> coefficients
  Fneu = params.conv_flux_ptr([],[],params);
  FNneu = Fneu; % simple identical coefficients!
end;
%| \docupdate 
