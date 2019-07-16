function Umean = fv_element_mean(model, model_data, U,I)
% function Umean = fv_element_mean(model, model_data, U,I)
% function computing the element averages of a discrete function U
% in the grid elements with indices I. Most arguments are dummy,
% but required for more general discrete functions, e.g. p1, etc.

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


% Bernard Haasdonk 16.5.2007

Umean = U(I);

