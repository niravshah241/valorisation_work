function reduced_data = gen_reduced_data(model, detailed_data)
%function reduced_data = gen_reduced_data(model, detailed_data)
%
% function computing reduced_data from model and detailed_data 
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

reduced_data = model.gen_reduced_data(model, detailed_data);

