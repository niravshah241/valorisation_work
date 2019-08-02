function detailed_data = gen_detailed_data(model, model_data)
%function detailed_data = gen_detailed_data(model, model_data)
%
% function computing detailed_data from model and model_data
% simple evaluation of function pointer in model

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

detailed_data = model.gen_detailed_data(model, model_data);
 