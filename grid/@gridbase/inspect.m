function inspect(grid, params)
%function inspect(grid, params)
% function plotting various properties of the current grid.
%
% The following plots are generated:
% - element numbers
% - the neighbour types and
% - neighbour indices
% - cog-coordinates
% - edge-midpoint-coordinates
%
% Parameters:
%  params: structure holding control fields for the functions
%          plot_element_data() and plot_sequence() used by this method.

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

%set element-data plot function for sequence plots

params.title = 'element indices';
params.plot = @plot_element_data;
plot_sequence((1:grid.nelements)',grid,params);
%plot_fv_data([1:grid.nelements; 1:grid.nelements]', params, ...
%	     'element indices');

params.title = 'boundary types';
plot_sequence(grid.NBI.*(grid.NBI<0),grid,params);
%plot_fv_data(grid.NBI.*(grid.NBI<0), params, ...
%	     'boundary types');

params.title = 'neighbour indices';
plot_sequence(grid.NBI,grid,params);
%plot_fv_data(grid.NBI, params, ...
%	     'neigbour indices');

params.title = 'cog coordinates';
plot_sequence([grid.CX,grid.CY],grid,params);
%plot_fv_data([grid.CX,grid.CY], params, ...
%	     'cog coordinates');

params.title = 'edge midpoint x-coordinates';
plot_sequence(grid.ECX,grid,params);
%plot_fv_data([grid.ECX], params, ...
%	     'edge midpoint x-coordinates');

params.title = 'edge midpoint y-coordinates';
plot_sequence(grid.ECY,grid,params);
%plot_fv_data([grid.ECY], params, ...
%	     'edge midpoint y-coordinates');

% additional check of consistency
%res = check_consistency(grid);

figure;
subplot(2,2,1),plot(grid),hold on, plot(grid.CX,grid.CY,'or');
title('centroids')

subplot(2,2,2),plot(grid),hold on, plot(grid.SX,grid.SY,'ob');
title('circumcenters')

subplot(2,2,3),plot(grid),hold on, plot(grid.ECX(:),grid.ECY(:),'og');
title('edge centroids')

subplot(2,2,4),plot(grid),hold on, plot(grid.ESX(:),grid.ESY(:),'ok');
title('circumcenter-edge intersections')


