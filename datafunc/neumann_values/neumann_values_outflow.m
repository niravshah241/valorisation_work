function FNneu = neumann_values_outflow(glob,U,normals,params)
%function FNneu = neumann_values_outflow(glob,U,normals,params)
%
% function computing neuman-values by pointwise evaluation. 
%
% glob: columnwise coordinate vectors of global points to be evaluated
% Uneu: columnwise U values in points in case of a u-dependent flux
% normals: Columnwise corresponding unit normal vectors
% 
% FNneu = F(Uneu, Xneu, Yneu) * n
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

Fneu = params.conv_flux_ptr(glob,U,params);
if decomp_mode == 0
  FNneu = Fneu(:,1) .* normals(:,1) + Fneu(:,2) .* normals(:,2);
elseif decomp_mode == 1
  Q_Fneu = length(Fneu);
  FNneu = cell(Q_Fneu,1);
  FNneu{:} = zeros(size(glob,1),1);
  for q = 1:Q_Fneu
    FNneu{q} =  Fneu{q}(:,1).*normals(:,1) + Fneu{q}(:,2).*normals(:,2);
  end;
else % mode 2
     % error('Not sure what to do in this case.')
     FNneu = Fneu;  % perhaps this is right?
end;
%| \docupdate 
