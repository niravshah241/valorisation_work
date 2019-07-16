function [flux, lambda] = conv_flux_velocity_matrixfile(glob,U,params)
%function [flux[, lambda]] = conv_flux_velocity_matrixfile(glob,U,params)
% function computing the convective flux `f(u)` of a convection problem.
%
% required fields of params:
% t: real time value in case of time-dependent flux  
% lambda : factor for velocity computation: `v = - \lambda \nabla p` a velocity
%          field is used, which is computed from an elliptic problem for the
%          pressure by compute_pressure_gdl2() (some gdl-elements in a row,
%          Neumann-0 and dirichlet conditions between 4 and 1 for the pressure,
%          linear decreasing along the channel) (reasonable domain
%          '[0,1e-3] x [0,0.25e-3]')
% velocity_matrixfile: filename of the mat-file storing the velocity field
%
% optional fields of params:
% filecache_velocity_matrixfile_extract: a boolean value indicating whether a
% cached version of the velocity field exists.
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

if decomp_mode == 2 % i.e. 'coefficients'
  flux   = 1; % factor one for single component
  lambda = [];
elseif decomp_mode < 2 % i.e. 'none' or 'components'
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% loading flux matrix.
  % load V.X: list of X-Koordinaten
  % load V.Y: list of Y-Koordinaten
  % load V.Vy: list of Y-component of velocity 
  % load V.Vx: list of Y-component of velocity 
  % load V.lambda: lambda for CFL computation 
  
  fullfn = fullfile(rbmatlabhome,'datafunc','mat',...
    params.velocity_matrixfile); 
  
  if ~exist(fullfn,'file') && ~cache('exist',fullfn)
    error(['warning: velocity_matrixfile not existing. ',...
      'Call gen_velocity_matrixfile with suitable' ...
      ' divclean_mode.']);
  end;
   
  % the following file extraction can be expensive, so the
  % extraction can additionally be cached.

  if ~params.filecache_velocity_matrixfile_extract % expensive call
    [flux_lin.Vx, flux_lin.Vy, lambda] = ...
      velocity_matrixfile_extract(fullfn,glob(:,1),glob(:,2));
  elseif params.filecache_velocity_matrixfile_extract==1 
    % cached functioncall
    [flux_lin.Vx, flux_lin.Vy, lambda] = ...
      filecache_function(@velocity_matrixfile_extract, ...
      fullfn,glob(:,1),glob(:,2));
  else % mode 2, i.e assume, that velocityfilename is set correctly to
    % requested points
    %	   V = load_cached(fullfn);
    V = cache('load',fullfn);
    % find integer indices such that V.X(j)==X(i) and V.Y(j) = Y(i) 
    % if whole flux-matrix is requested
    if ~isequal(glob(:,1),V.X(:)) || ~isequal(glob(:,2),V.Y(:))
      error('matrixfile does not fit to requested points!');
    end;
    % i = 1:length(X);
    flux_lin.Vx = V.Vx;
    flux_lin.Vy = V.Vy;
    %NaN*ones(size(X));
    %	   flux_lin.Vy = NaN*ones(size(X));
    %	   flux_lin.Vx(i) = V.Vx(j);
    %	   flux_lin.Vy(i) = V.Vy(j);
    lambda = V.lambda;
  end;

  flux = [ flux_lin.Vx.*U(:), flux_lin.Vy.*U(:) ];

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end of loading flux matrix.
  if decomp_mode == 1 % put single matrix into cell array
    flux = {flux};
  end;
else
  error('unknown decomp_mode');
end

