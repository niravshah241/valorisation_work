function [X_trans, Y_trans] = geometry_transformation(X,Y,params)
%function [X_trans,Y_trans] = geometry_transformation(X,Y,params)
%
% function applying a geometry transformation function on the
% coordinates X,Y. The type of geometry manipulation is given through
% the params structure.
%
% required fields of params:
%  geometry_transformation  : 'spline', 'none'
%  geometry_transformation_spline_x,
%  geometry_transformation_spline_y : coordinates of interpolation for
%                                     spline line in case of 'spline'
% 
% optional fields of params:
%  hill_height              : if set it overwrites the y_hill(2)
%
% The different transformations:
%  'spline' transforms the upper line into a spline given by
%  interpolation points. The points beneath this "hill" are linearly
%  interpolated.
%
% Martin Drohmann 22.02.2008

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


if isstruct(params) && ~isfield(params,'geometry_transformation')
  params.geometry_transformation = 'none';
end

if isequal(params.geometry_transformation, 'none')
  X_trans = X;
  Y_trans = Y;
elseif isequal(params.geometry_transformation, 'spline')
  grid_y_min  = params.yrange(1);
  grid_height = params.yrange(2) - params.yrange(1);

  cs = spline_select(params);
  
%  transfunc_Y = @(X,Y) ( 1/grid_height + ppval(cs, X) ) ...
%                        .* (Y - grid_y_min) ...
%               		  + grid_y_min;
  transfunc_Y = @(X,Y) (Y-grid_y_min) ./ grid_height ...
                       .* (grid_height + ppval(cs,X)) ...
                       + grid_y_min;
  X_trans = X;
  Y_trans = transfunc_Y(X,Y);
else
    error(['unknown method for geometry transformation: ', params.geometry_transformation]);
end

%| \docupdate 
