function str = disp_deprecated(varargin)
%function str = disp_deprecated([alternative])
%
% determine name of caller and print deprecated message plus
% eventual alternative

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


% Bernard Haasdonk 4.9.2009

s = dbstack;
fname = s(2).name;
disp(['warning: ',fname,' deprecated!!!']);
if (length(varargin)>0)
  disp(['         ',varargin{1},' should be used instead.']);
end;
%| \docupdate 
