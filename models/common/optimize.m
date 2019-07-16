function [opt_data,model] = optimize(model, model_data, detailed_data, reduced_data)
%opt_data = optimize(model, model_data, detailed_data, reduced_data)
%
%Method starting optimization and returning results in opt_data.
%
% Markus Dihlmann 28.04.2010
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


[opt_data,model] = model.optimization.optimizer(model, model_data, detailed_data, reduced_data);

