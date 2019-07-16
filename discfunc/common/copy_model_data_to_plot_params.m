function plot_params = copy_model_data_to_plot_params(model, plot_params)
% function plot_params = copy_model_data_to_plot_params(model, plot_params)
% Helper function copying extracting relevant information for plot_params from
% the model.
%
% This is especially for use with models using geometry transformation. Here,
% the geometry parameters also need to be known by the 'plot_params'
%
% Parameters:
%  plot_params: structure holding the parameters for plotting.
%
% Return values:
%  plot_params: the updated structure enriched with fields extracted from the
%  'model'.

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


plot_params.geometry_transformation = model.geometry_transformation;

plot_params.xrange                           = model.xrange;
plot_params.yrange                           = model.yrange;
if isfield(model, 'geometry_spline_type')
  plot_params.geometry_spline_type             = model.geometry_spline_type;
  plot_params.geometry_transformation_spline_x = model.geometry_transformation_spline_x;
  plot_params.geometry_transformation_spline_y = model.geometry_transformation_spline_y;
  plot_params.hill_height                      = model.hill_height;
end

if isfield(model, 'postprocess')
  plot_params.postprocess                    = model.postprocess;
  plot_params.gravity                        = model.gravity;
  plot_params.clim                           = model.clim;
end

