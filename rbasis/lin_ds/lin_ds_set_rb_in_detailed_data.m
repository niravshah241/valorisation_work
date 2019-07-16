function detailed_data = ...
    lin_ds_set_rb_in_detailed_data(detailed_data,RB)
%function detailed_data = lin_ds_set_RB_in_detailed_data(detailed_data,RB)
%
% function setting fields V=W=RB in detailed_data

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

% Bernard Haasdonk 21.9.2009

detailed_data.V = RB;
if ~isempty(RB) 
    detailed_data.W = detailed_data.G * RB;
else
    detailed_data.W = [];
end;
