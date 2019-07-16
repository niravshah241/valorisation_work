function s = copy_field_if_exists(s,source,fieldname,defaultvalue)
%function s = copy_field_if_exists(s,source,fieldname,defaultvalue)
%
% function setting the field fieldname in the struct s. Either the
% field is copied from source if it exists, or the defaultvalue is set

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


% Bernard Haasdonk 19.3.2010

if isfield(source,fieldname)
  s = setfield(s,fieldname,getfield(source,fieldname));
else
  s = setfield(s,fieldname,defaultvalue);
end;