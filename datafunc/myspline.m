function [p_mu] = myspline(x_hill, y_hill)

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


p_mu = spline(x_hill, y_hill);
%p_mu = mkpp([0 1], [-y_hill(2) y_hill(2)]);

% TO BE ADJUSTED TO NEW SYNTAX
%| \docupdate 
