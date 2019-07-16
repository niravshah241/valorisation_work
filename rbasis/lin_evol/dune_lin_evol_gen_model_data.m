function model_data = dune_lin_evol_gen_model_data(model)
%function model_data = dune_lin_evol_gen_model_data(model)
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


model.mexptr('gen_model_data');

% model data is saved in mex process, the following is just for compatibility
% with RBmatlab models.
model_data        = [];
model_data.W      = [];
model_data.grid   = 'initialized in mexfile';
model_data.mexptr = model.mexptr;

%| \docupdate 
