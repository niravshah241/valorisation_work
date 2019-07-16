function [scalar_dofs, scalar_df_info] = fem_scalar_component(dofs,ncomp,df_info)
%function [scalar_dofs, scalar_df_info] = fem_scalar_component(dofs,ncomp,df_info)
%
% extract single component ncomp of vectorial discrete function and 
% generate new scalar fem function dof vector of same degree.
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


% Bernard Haasdonk 12.1.2011

indices = (0:(df_info.ndofs/df_info.dimrange-1))*df_info.dimrange + ncomp;
scalar_dofs = dofs(indices)';
scalar_params = [];
scalar_params.dimrange = 1;
scalar_params.pdeg = df_info.pdeg;
% perhaps need more in scalar_params
scalar_df_info = feminfo(scalar_params,df_info.grid);