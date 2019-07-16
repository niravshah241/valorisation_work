function res = isequal_ignore_handles(s1,s2)
%function res = isequal_ignore_handles(s1,s2)
%
% generalization of isequal for structures with possibly
% function_handles. as isequal does not work for these,
% function_handles are skipped
%
% parameters:
%  s1: general variable to be compared
%  s2: general variable to be compared
%
% return value:
%  boolean, indicating wether the two variables are identical (while ignoring
%  function handles)

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


res = 0;

if ~isequal(class(s1),class(s2))
  return;
end;
  
switch class(s1)
 case 'function_handle'
  res = 1; 
  return
 case {'double','char','single','logical','rectgrid','triagrid','onedgrid'}
  res = isequal(s1,s2);
 case 'cell'
  if length(s1)~=length(s2)
    return;
  else
    res = 1;
    for i=1:length(s1)
      res = res && isequal_ignore_handles(s1{i},s2{i});
    end;
  end;
 case 'struct'
  fn1 = fieldnames(s1);
  fn2 = fieldnames(s2);
  if length(intersect(fn1,fn2))~=length(fn1)
    res = 0;
    return;
  end;
  res = 1;
  for i = 1:length(fn1)
    res = res && ...
	  isequal_ignore_handles(s1.(fn1{i}),s2.(fn1{i}));
  end;
 otherwise
  error('class type not yet supported');
end;

