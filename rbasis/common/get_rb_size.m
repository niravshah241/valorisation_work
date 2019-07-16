function rb_size = get_rb_size(model,detailed_data)
%function rb_size = get_rb_size(model,detailed_data)
%
% Function returning the size of the actual basis in detailed_data
%
% Required fields of detailed_data:
%   RB: The reduced basis
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


if isfield(detailed_data,'RB')
  rb_size = size(detailed_data.RB,2);
else
  rb_size = 0;
end


