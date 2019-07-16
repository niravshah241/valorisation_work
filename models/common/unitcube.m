function model=unitcube(model)
% function model=unitcube(model)
% function adding fields to model for generating a 2D rectgrid with '100 x 100'
% elements on the unit-square

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


model.gridtype      = 'rectgrid';
model.xrange        = [0,1];
model.yrange        = [0,1];
model.xnumintervals = 100;
model.ynumintervals = 100;

