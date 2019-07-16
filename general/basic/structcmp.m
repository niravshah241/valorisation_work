function [iseq, diffs, add1, add2] = structcmp(s1,s2,ignorelist)
%function [iseq, diffs, add1, add2] = structcmp(s1,s2,ignorelist);
%
% compare two structures s1 and s2 and ignoring fields given in the
% cell-array of fielnames ignorelist.
% iseq = 1 if the structures are identical except ignorelist
% diffs is a cell array of the fieldnames persistent in both
% structs but differing, add1 are the fields existing in s1 but
% not in s2, add2 is a cell array of fieldnames existing in
% s2 but not in s1
% function handles are ignored automatically, as they cannot be compared

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


% Bernard Haasdonk 23.5.2007

diffs = {};

if (nargin<3) || isempty(ignorelist)
  ignorelist = {};
end;

f1 = fieldnames(s1);
f1 = setdiff(f1, ignorelist);
f2 = fieldnames(s2);
f2 = setdiff(f2, ignorelist);
add1 = setdiff(f1,f2);
add2 = setdiff(f2,f1);
fi = intersect(f1,f2);

for i = 1:length(fi)
  if isa(s1.(fi{i}), 'function_handle')
    if ~isequal(func2str(s1.(fi{i})), func2str(s2.(fi{i})))
      diffs = [ diffs; {fi(i)}];
    end
  elseif iscell(s1.(fi{i})) ...
           && all(cellfun(@(x) isa(x, 'function_handle'), s1.(fi{i})))
    strings1 = cellfun(@func2str, s1.(f1{i}), 'UniformOutput', false);
    strings2 = cellfun(@func2str, s1.(f1{i}), 'UniformOutput', false);
    if ~isequal(strings1, strings2)
      diffs = [ diffs; {fi(i)} ];   
    end
  elseif ~isequal(s1.(fi{i}),s2.(fi{i}))
    diffs = [diffs; {fi(i)}];
  end
end;

iseq = 1;
if ~isempty(diffs) ||  ...
      ~isempty(add1) || ...
      ~isempty(add2)  
  iseq = 0;
end;
%| \docupdate 
