function rb_sim_data = rb_simulation(model,reduced_data)
%function rb_sim_data = rb_simulation(model,reduced_data)
%
% method performing reduced simulation.
% simple call of pointer in model

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


% Martin Drohmann and Bernard Haasdonk 21.7.2009

rb_sim_data = model.rb_simulation(model, reduced_data);

