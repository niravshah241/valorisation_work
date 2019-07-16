function rb_sim_data = lin_ds_rb_reconstruction(model,detailed_data, rb_sim_data)
%function rb_sim_data = lin_ds_rb_reconstruction(model,detailed_data, rb_sim_data)
% 
% function computing a reconstruction of the reduced simulation
% trajectory, i.e. simply multiplication of reduced coefficients with reduced basis

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


% Bernard Haasdonk 2.4.2009

N = size(rb_sim_data.Xr,1);
rb_sim_data.X = detailed_data.V(:,1:N) * rb_sim_data.Xr;
