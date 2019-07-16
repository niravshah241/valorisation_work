function p = plot_grid(grid,params)
%function p = plot_grid(grid,params)
%
% plot method for a cubegrid. A line plot is performed.
% Text can be plotted in the leafs of the grid.
%
% p is the list of handles to the graphics primitives
%
% optional params:
%   color : RGB vector
%   shrink_factor: if this flag is given, the elements are plotted shrinked
%   plot_level:   if this flag is nonzero, a level plot is performed
%   level  : this integer must be specified in case of level-plot    
%   plot_patch    : if this flag is set (only 2D and 3D) the plot is done
%                   by colored patches
%   axis_equal    : if this flag is set, set axis to equal scale
%
%   leaftext      : vector of strings. Length must be equal to number of
%                   leafelements. The text is shown in the cogs (center of
%                   gravities) of the leafelements.
%
% Markus Dihlmann 18.02.2010

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


p=plot(grid,params);

cogs = get(grid,'leafcogs');

if isfield(params,'text')
if(length(params.text)==get(grid,'nleafelements'))
    for i=1:length(cogs(:,1))
    text(cogs(i,1)-0.025,cogs(i,2),params.text{i},'HorizontalAlignment','left','FontSize',14,'FontWeight','b');
    end
else
    disp('Couldnt show text in plot because text vector has not the appropriate length regarding nleaf_elemnts');
end  
end
end
