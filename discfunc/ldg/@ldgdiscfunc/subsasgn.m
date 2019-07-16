function res = subsasgn(df, S, val) 
%function res = subsasgn(df, S, val) 
%
% function allowing access to member-variables of df
% this simplifies use remarkably, although also some 
% private member-variables can be read, so data-capsulation is 
% slighlty violated by this. lternative would be to provide around 
% 20 methods for access to all data fields. Therefore this is the
% better compromise.

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


% Bernard Haasdonk 9.5.2007

res = builtin('subsasgn', df, S, val);%| \docupdate
