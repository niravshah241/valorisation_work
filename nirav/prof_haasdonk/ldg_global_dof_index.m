function gids = ldg_global_dof_index(params, grid)
%function gids = ldg_global_dof_index(params, grid)
%
% function computing the global dof indices of a ldg function space
% params can either be a df function or a df_info or a simple
% struct with field ndofs, ndofs_per_element
% gids is a nelements x nlocaldofs matrix 
% where gids(elid, ldofid) gives the global dof index of the local
% basis function number ldofid on element with number elid in the
% current grid.

% Bernard Haasdonk 12.8.2017
gids = reshape(1:params.ndofs, params.ndofs_per_element, grid.nelements)';

