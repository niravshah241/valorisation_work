function homestr = rbmatlabhome
%  function homestr = rbmatlabhome
%
% function returning the home environment variable pointing to the 
% RBmatlab-subdirectory 

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

  
% Bernard Haasdonk 21.7.2006
 
 %%%  the following call is very inefficient, i.e. 1400 calls = 4 sec!!!! 
 %homestr = [fileparts( which('startup_rbmatlab')),filesep]; 

 %%% better: environment-variable method: fraction of a second for
 %1000 calls. good.
 
 homestr = getenv('RBMATLABHOME');
 
 %  homestr = getenv('HOME');
 %  if isequal(homestr(1:3),'C:\')
 %    homestr = [homestr,'sync_lcars'];
 %  end;
 %  homestr = fullfile(homestr,'matlab','RBmatlab');
  
% TO BE ADJUSTED TO NEW SYNTAX
%| \docupdate 
