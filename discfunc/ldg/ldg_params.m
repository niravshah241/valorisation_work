function params = ldg_params(nelements,pdeg,dimrange,params)
% function setting default parameter values in params.
%
% params is extended by corresponding fields.

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


% Bernard Haasdonk 18.1.2010

if nargin < 4
  params = [];
end;
params.nelements = nelements; 
params.pdeg = pdeg;
params.dimrange = dimrange; % vectorial function
params.ndofs = ldg_ndofs(params);
params.ndofs_per_element = ldg_ndofs_per_element(params);
%| \docupdate 
