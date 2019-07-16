function p = rb_plot_output_estimation(simulation_data,model)
%function p = rb_plot_output_estimation(simulation_data,model)
%
% function plotting the estimated output quantity and its error bounds
%
% required fileds of simulation_data:
%   Delta_s: row vector of output error data
%

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


% Bernard Haasdonk 16.5.2008


s = simulation_data.s;
s_up  = s + simulation_data.Delta_s;
s_low  = s - simulation_data.Delta_s;
plot([s,s_up,s,s_low]);
xlabel('time step');
ylabel('output + bounds');
title(['output estimation']);
