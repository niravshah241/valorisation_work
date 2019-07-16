function [p_mu] = spline_select(model)
%function [p_mu] = spline_select(model)
%
% This function is a wrapper for spline functions that are used for geometrical
% transformations of the grid.

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


if isequal(model.geometry_spline_type, 'cubic')

  x_hill  = model.geometry_transformation_spline_x;
  y_hill  = model.geometry_transformation_spline_y;
  if model.hill_height >= 0
    y_hill(2) = model.hill_height;
  end

  p_mu = spline(x_hill, y_hill);

elseif isequal(model.geometry_spline_type, 'affine')

%  x_hill  = model.geometry_transformation_spline_x;
  y_hill  = model.geometry_transformation_spline_y;
  if model.hill_height >= 0
    y_hill(2) = model.hill_height;
  end

%  p_mu = mkpp(x_hill, [(y_hill(2)-y_hill(1))/(x_hill(2)-x_hill(1)),y_hill(1);
%                       (y_hill(3)-y_hill(2))/(x_hill(3)-x_hill(2)),y_hill(2)]);
   p_mu = mkpp([0 1], [-y_hill(2) y_hill(2)]);
%  plot(ppval(p_mu, linspace(0,1)));
end

%| \docupdate 
