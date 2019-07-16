function p = lin_ds_plot_sim_data_output(model,model_data,sim_data,params)
%function p = lin_ds_plot_sim_data_output(model,model_data,sim_data[,params])
%
% plot of dynamical system simulation output results
%
% Return values:
%  - p: plot handles

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


% Bernard Haasdonk 7.9.2009

if ~isfield(params,'output_plot_indices')
  output_plot_indices = 1:model.dim_y;
else
  output_plot_indices = params.output_plot_indices;
end;

if ~isfield(params,'plot_output_title')
  plot_output_title = 'output';
else
  plot_output_title = params.plot_output_title;
end;

figure, 
p = plot(sim_data.time,sim_data.Y(output_plot_indices,:));
legend(num2str((output_plot_indices(:))));
title(plot_output_title);
xlabel('time t');
ylabel('y_i(t)');
axis tight;
set(p,'Linewidth',2)
 
