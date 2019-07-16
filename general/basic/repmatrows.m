function rmat = repmatrows(A,times)
%
% function stretching a matrix in first direction by 
% copying the rows

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


% Bernard Haasdonk 19.7.2006

rmat = reshape(repmat(A',times,1),...
	       size(A,2),size(A,1)* times)';
% TO BE ADJUSTED TO NEW SYNTAX
%| \docupdate 
