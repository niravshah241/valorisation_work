function fname = cache_velocity_matrixfile_extract(model,...
                                                   X, Y, ext)
%function fname = cache_velocity_matrixfile_extract(model, ...
%                                                  X, Y, ext)
%
% function checking the availability of the velocity matrixfile in the cache
% required for evaluation of the velocity in the given points. 
% No real files are generated, but the persitent cache variable is used
% for storing these files.

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

% Bernard Haasdonk 4.9.2007

%fname = fullfile(rbmatlabhome,'datafunc','mat',[ext,'_', ...
%		    params.velocity_matrixfile]);

fname = [ext,'_', model.velocity_matrixfile];
fullfname = fullfile(rbmatlabhome,'datafunc','mat',fname);

%if ~exist(fullfname,'file')
if ~cache('exist',fullfname)
  if model.verbose > 9
    disp(['generating velocity data in cache: ',fname]);
  end;
  completefn = fullfile(rbmatlabhome,'datafunc','mat',...
			model.velocity_matrixfile);
  [Vx,Vy, lambda] = ...
      velocity_matrixfile_extract(completefn,X, Y);
  tmp.X = X;
  tmp.Y = Y;
  tmp.lambda = lambda; 
  tmp.Vx = Vx;
  tmp.Vy = Vy;  
  %  save(fullfname,'X','Y','lambda','Vx','Vy');
  cache('save',fullfname,tmp);
end;

%| \docupdate 
