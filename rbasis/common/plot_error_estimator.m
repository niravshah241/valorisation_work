function p = plot_error_estimator(offline_data,params)
% function p = plot_error_estimator(offline_data,params)
%
% function plotting the error-estimator landscape over the
% parameter space.
%
% required fields of params:
%
%     mu_ranges: cell array of the mu_ranges
%     num_plot_intervals : vector indicating the number of
%                  intervals for the rb_simulation test. If not
%                  set, default [10],[10 10],[10 10 10] is taken
%                  depending on the dimensionality of the parameter space

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


% Bernard Haasdonk 6.6.2007

if ~isfield(params,'num_plot_intervals')
  params.num_plot_intervals = 10 * ones(length(params.mu_names),1);
end;

par.range = params.mu_ranges;
par.numintervals = params.num_plot_intervals;
plotgrid = cubegrid(par);

M = get(plotgrid,'vertex')';

params.error_indicator = 'estimator';
error_estimators = rb_test_indicator([],offline_data,M,[],params);

plot_leafvertex_data(error_estimators,plotgrid,params);

% TO BE ADJUSTED TO NEW SYNTAX
%| \docupdate 
