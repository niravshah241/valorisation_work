function [flux, lambda] = conv_flux_gdl2(glob,U,params)
%function [flux[, lambda]] = conv_flux_gdl2(glob,U,params)
% Function computing the convective flux `f(u)` of a convection problem.
%
% required fields of params:
% t: real time value in case of time-dependent flux  
% lambda : factor for velocity computation: `v = - \lambda \nabla p`
%            a velocity field is used, which is computed from an
%            elliptic problem for the pressure by compute_pressure_gdl2() (some
%            gdl-elements in a row, Neumann-0 and dirichlet conditions between
%            4 and 1 for the pressure, linear decreasing along the channel)
%            (reasonable domain '[0,1e-3] x [0,0.25e-3]')
%

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

% glob column check
if params.debug
  if ~isempty(glob) && size(glob,1) < size(glob,2)
    warning('coordinates in variable glob are given row-wise, but expected them to be column-wise');
    if params.debug > 2
      keyboard;
    end
  end
end

% determine affine_decomposition_mode as integer  
decomp_mode = params.decomp_mode;


if isempty(glob)
  X = [];
  Y = [];
else
  X = glob(:,1);
  Y = glob(:,2);
end

% load pressure from file and generate velocity field
p = load(fullfile(rbmatlabhome,'datafunc','mat','gdl_pressure2.mat'),...
  'p','params');
[ux, uy] = evaluate_gradient(X,Y,p.p,p.params);  
%if iscell(ux)
%  lin_flux = cellfun(@(X,Y)([-X*params.lambda, -Y*params.lambda]), ...
%                             ux, uy, 'UniformOutput', false);
%else
  lin_flux = [-ux*params.lambda, -Y*params.lambda];
%end

if decomp_mode == 2
  flux = 1;
  lambda = [];
elseif decomp_mode < 2
  % inverse of maximum absolute velocity value 
  % correct up to a constant factor
  lambda = max(max(lin_flux))^-1; 
  if ~isempty(lin_flux) % otherwise size 0/0 => 0/1 !!!!!
    flux = [ lin_flux(:,1).* U(:), lin_flux(:,2).* U(:) ];
  else
    flux = ones(0,2);
  end;
  if decomp_mode == 1
    flux = { flux };
  end
else
  error('unknown decomp_mode');
end


