function model = set_mu_in_model_and_base_model_old(model,mu)
%function model = set_mu_in_model_and_base_model_old(model,mu)
% 
% function setting the parameter vector mu in the model and the
% submodel base_model

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

  
% Bernard Haasdonk 12.2.2010

model = model.set_mu(model,mu);
model.base_model = model.base_model.set_mu(model.base_model,mu);


