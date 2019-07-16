function p = rb_plot_reconstruction(detailed_data,simulation_data,model)
%function p = rb_plot_reconstruction(detailed_data,simulation_data,model)
% 
% (trivial) function computing a detailed reconstruction by linear
% combination of the coefficients in the simulation data with the
% orthonormal reduced basis RB and plotting

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

model.title = 'RB-solution reconstruction';
Uapprox = detailed_data.RB(:,1:size(simulation_data.a,1)) * ...
	  simulation_data.a;
plot_element_data_sequence(detailed_data.grid, Uapprox, model);
