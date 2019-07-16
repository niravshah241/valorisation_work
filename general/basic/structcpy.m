function [s1] = structcpy(s1,s2)
% function [s1] = structcpy(s1,s2)
%
% copies the fields of structure s2 into structure s1. If the field to be
% copied does not exist in s1 yet, a field with the appropriate name is
% created.

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


% Martin Drohmann, 01.09.2009

f2 = fieldnames(s2);
for i = 1:length(f2)
  s1.(f2{i}) = s2.(f2{i});
end

%| \docupdate 
