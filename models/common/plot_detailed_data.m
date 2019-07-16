function plot_detailed_data(model, detailed_data, plot_params)
%function plot_detailed_data(model, detailed_data, plot_params)
%
% function plotting the detailed data in a problem specific manner
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


% Bernard Haasdonk 26.8.2009
%

if nargin <= 2
  plot_params = [];
end

model.plot_detailed_data(model,detailed_data,plot_params);


