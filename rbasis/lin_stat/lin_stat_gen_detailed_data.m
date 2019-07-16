function detailed_data = lin_stat_gen_detailed_data(model,model_data)
%function detailed_data = lin_stat_gen_detailed_data(model,model_data)
%
% function computing the high-dimensional parameter independent 
% data required for rb-treatment of a fem-problem, comprises all
% expensive computations

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


% B. Haasdonk 22.2.2011

detailed_data = [];
detailed_data.grid = model_data.grid;
detailed_data.df_info = model_data.df_info;
detailed_data = rb_basis_generation(model,detailed_data);
