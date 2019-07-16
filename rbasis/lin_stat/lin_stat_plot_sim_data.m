function p = lin_stat_plot_sim_data(model,model_data,sim_data,plot_params);
%function p = lin_stat_plot_sim_data(model,model_data,sim_data,plot_params);

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


% B. Haasdonk 22.2.2011

% plot result
if nargin < 4
  plot_params = [];
end;
if ~isfield(plot_params,'subsampling_level');
  plot_params.subsampling_level = 10;
end;
if ~isfield(plot_params,'title');
  plot_params.title = '';
end;
p = plot(sim_data.uh,plot_params);
title(plot_params.title);
% later: output as text label or as screen print, etc.

  