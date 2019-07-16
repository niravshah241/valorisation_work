function reduced_data_subset = extract_reduced_data_subset(model, reduced_data)
%function reduced_data_subset = extract_reduced_data_subset(model, reduced_data)
%
% method which modifies reduced_data, which is the data, that will
% be passed to the online-simulation algorithm. 
% Typically, this routine only does a submatrix extraction of the Mmax, Nmax 
% sized offline-objects to produce M, N sized objects for the real simulation.
% Required fields of model:
%  N       : number of reduced basis vectors to choose

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


reduced_data_subset = model.reduced_data_subset(model, reduced_data);
%| \docupdate 
