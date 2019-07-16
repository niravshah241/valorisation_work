function [scalar_dofs, scalar_params] = ldg_scalar_component(dofs,ncomp,params)
%function [scalar_dofs, scalar_params] = ldg_scalar_component(dofs,ncomp,params)
%
% extract single component ncomp of vectorial discrete function and 
% generate new scalar ldg function dof vector of same degree.
%
% params must provide params.pdeg, params.nelements, params.dimrange
% params.ndofs

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


% Bernard Haasdonk 27.1.2009

indices = (0:(params.ndofs/params.dimrange-1))*params.dimrange + ncomp;
scalar_dofs = dofs(indices)';
scalar_params = [];
scalar_params.dimrange = 1;
scalar_params.pdeg = params.pdeg;
scalar_params.nelements = params.nelements;
scalar_params.ndofs_per_element = ldg_ndofs_per_element(scalar_params);
scalar_params.ndofs = ldg_ndofs(scalar_params);
%| \docupdate 
