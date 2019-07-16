function Mstr = matrix2str(M)
%function Mstr = matrix2str(M)
%
% function generating a string (matlab format) of a given double
% matrix M

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


% Bernard Haasdonk 30.1.2009

fstr = [repmat('%30.18d, ',1,size(M,2)-1),'%30.18d; ...\n'];
Mstr = ['[ ',sprintf(fstr,M'), '];'];%| \docupdate 
