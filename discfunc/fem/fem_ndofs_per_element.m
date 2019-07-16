function res = fem_ndofs_per_element(params)
%function res = fem_ndofs_per_element(params)
%
% function computing number of dofs per element based on 
% params.dimrange and params.pdeg.

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


% Bernard Haasdonk 11.1.2011

res = params.dimrange*(params.pdeg+1)*(params.pdeg+2)/2;
 %| \docupdate 
