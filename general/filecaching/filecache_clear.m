function filecache_clear
%function filecache_clear
%
% function clearing the filecache. should be called at each new
% program start which uses filecaching. But at least of course before
% each new time-measurement!!!

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


% Bernard Haasdonk 22.5.2007

d = dir(filecache_path);
for i = 1:length(d)
  if ~ismember(d(i).name,{'.','..'})
    delete(fullfile(filecache_path,d(i).name));
  end;
end;


% TO BE ADJUSTED TO NEW SYNTAX
%| \docupdate 
