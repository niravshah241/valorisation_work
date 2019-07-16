function merged_plot_params = merge_model_plot_params(model, plot_params)
%function merged_plot_params = merge_model_plot_params(model, plot_params)
%
% function setting some default values in plot_params and copying
% fields from model, then from plot_params, hence plot_params 
% can overload settings in the model.

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


% B. Haasdonk 11.6.2010

fieldnames = {'plot',...
	      'plot_slice',...
	      'geometry_transformation',...
	      'xrange',...
	      'yrange',...
	      'geometry_spline_type',...
	      'geometry_transformation_spline_x',...
	      'geometry_transformation_spline_y',...
	      'hill_height',...
	      'post_process',...
	      'gravity',...
	      'clim', ...
	      'range_lim', ...
	      'axis_equal', ...
	      'axis_tight', ...
	      'no_lines'...
	     };

merged_plot_params = plot_params;
for i = 1:length(fieldnames)
  fn = fieldnames{i};
  if isfield(model,fn)
    merged_plot_params = setfield(merged_plot_params,fn,getfield(model,fn));
  end;
  if isfield(plot_params,fn)
    merged_plot_params = setfield(merged_plot_params,fn,getfield(plot_params,fn));
  end;
end;

%| \docupdate 
