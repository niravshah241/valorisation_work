function sim_data = detailed_ei_simulation(model,model_data)
%function sim_data = detailed_ei_simulation(model,model_data)
%
% method performing detailed simulation with empirically interpolated operators
% models should have default parameters.
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

if nargin == 2
  params = [];
end

sim_data = model.detailed_ei_simulation(model, model_data);

%| \docupdate 
