function tmpstr = rbmatlabtemp
%  function tmpstr = rbmatlabtemp
%
% function returning the home environment variable at which the
% large space for data is available

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
  
tmpstr = getenv('RBMATLABTEMP');

%  homestr = getenv('HOME');
 %  if isequal(homestr(1:3),'C:\')
 %    homestr = [homestr,'sync_lcars'];
 %  end;
 %  homestr = fullfile(homestr,'matlab','RBmatlab');
  
% TO BE ADJUSTED TO NEW SYNTAX
%| \docupdate 
