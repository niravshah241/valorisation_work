function res = subsref(df, S) 
%function res = subsref(df, S) 
%
% function allowing access to member-variables of df
% this simplifies use remarkably, although also some 
% private member-variables can be read, so data-capsulation is 
% slighlty violated by this. lternative would be to provide around 
% 20 methods for access to all data fields. Therefore this is the
% better compromise.
%
% additionally, the operator() is overloaded: Now this is a local 
% evaluation of ldg functions. This enables identical syntax for
% evlauating analytical or discrete functions with arguments 
% (einds, lcoord, grid, params)

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

if S(1).type~='('
  res = builtin('subsref', df, S);
else
  res = evaluate(df,S.subs{:});
end;%| \docupdate 
