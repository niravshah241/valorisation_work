function model_data = init_model(model_name, varargin)
% function model_data = init_model(model_name, varargin)
%
% function that initializes the numerical name given by 'model_name'. This
% model needs to supply a M-File or a Mex-File named 'model_name'. For details
% see the README.modelinterface. It returns model specific data like
% information about the number and range of the parameters.

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


model_data = feval(model_name, 'init_model', varargin{:});

%| \docupdate 
