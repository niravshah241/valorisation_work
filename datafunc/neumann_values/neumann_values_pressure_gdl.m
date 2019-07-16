function FNneu = neumann_values_pressure_gdl(glob,U,normals,params)
%function FNneu = neumann_values_pressure_gdl(glob,U,normals,params)
%
% function computing neuman-values by pointwise evaluation. 
%
% glob: columnwise coordinate vectors of global points to be evaluated
% Uneu: columnwise U values in points in case of a u-dependent flux
% normals: Columnwise corresponding unit normal vectors
% 
% FNeu = no_flow on upper and lower boundary. left and
% right is a parabolic velocity profile with maximum params.c_neu_max.
% Uneu and Nxneu and Nyneu may be ignored.
% 'rightflow': outflow to right, noflow to upper and lower
%
% required fields of params:
%   params.c_neu_max : maximum flow value for 'pressure_gdl'
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

% set all to zero
FNneu = zeros(size(glob,1),1); 			
% left boundary negative velocity -> positive pressure gradient      
i = find(glob(:,1)<params.xrange(1)+eps);
FNneu(i) = params.c_neu_max * ...
    (params.yrange(2)-params.yrange(1))^(-2) * ...
    (glob(i,2)-params.yrange(1)).*(params.yrange(2)-glob(i,2))*4; 
% right boundary positive velocity -> negative pressure gradient    
i = find(glob(:,1)>params.xrange(2)-eps);
FNneu(i) = -params.c_neu_max * ...
    (params.yrange(2)-params.yrange(1))^(-2) * ...
    (glob(i,2)-params.yrange(1)).*(params.yrange(2)-glob(i,2))*4; 

if decomp_mode>0 
  error('function does not support affine decomposition!');
end;

%| \docupdate 
