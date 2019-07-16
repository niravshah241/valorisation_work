function rb_plot_detailed_data(model,detailed_data)
%function rb_plot_detailed_data(model,detailed_data)
% 
% function simply calling a problem specific implementation defined
% by the model.rb_problem_type, i.e.
% rb_plot_lin_evol_detailed_data, etc.

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


% Bernard Haasdonk 23.5.2007

plot_algorithm = ['rb_',params.rb_problem_type,...
		  '_plot_detailed_data'];

feval(plot_algorithm,detailed_data,params);

