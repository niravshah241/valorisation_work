function rb_sim_data = rb_reconstruction_default(model, detailed_data, rb_sim_data)
%function rb_sim_data = rb_reconstruction_default(model, detailed_data, rb_sim_data)
%
% (trivial) function computing a detailed reconstruction by linear
% combination of the coefficients in the simulation data with the
% orthonormal reduced basis RB

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


% Bernard Haasdonk 19.7.2006

rb_sim_data.U = detailed_data.RB(:,1:size(rb_sim_data.a,1)) * rb_sim_data.a;
