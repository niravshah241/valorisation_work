function res = cache(command, fn, varargin)
%function res = cache(command, fn, varargin)
%
% function simulating a ram-disk. 
%
% D = cache('load','myfile') :   loading of myfile from cache. If not
%                                in cache, the real file is loaded. Take
%                                care: if the file is changed while data is
%                                in cache, this is not detected!
% D = cache('clear','myfile') :  clearing myfile from cache. 
% D = cache('clear') :           the whole cache is cleared. 
% cache('save', 'myfile', data): saving data to myfile in
%                                cache, no real saving to disk is performed! 
% list = cache('list') :         return list of current entries
% res = cache('exist',fn) :      1 if dataset exists in cache, 0 otherwise
% cache('limit',limsize) :       sets the size limit to limsize (in byte)
%                                default: 1 MB

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

  
% Bernard Haasdonk 23.8.2007

persistent keys data limsize

if isempty(limsize)
  limsize = 10000000;
end;

%i = find(ismember(keys,fn));

switch command
 case 'exist'
  res = 0;
  
  for j=1:length(keys)
    if length(keys{j})==length(fn) & ...
	  isequal(keys{j},fn)
      res = 1;
      return;
    end;
  end;
  
 case 'load'

  i = [];
  for j=1:length(keys)
    if length(keys{j})==length(fn) & ...
	  isequal(keys{j},fn)
      i = j;
    end;
  end;
  
  if ~isempty(i) 
    res = data{i};
  else
    res = load(fn);
    data{end+1} = res;
    keys{end+1} = fn;
  end;
  
 case 'save'
  
  i = [];
  for j=1:length(keys)
    if length(keys{j})==length(fn) & ...
	  isequal(keys{j},fn)
      i = j;
    end;
  end;
  
  if ~isempty(i) 
    data{i} = varargin{1};
  else
    data{end+1} = varargin{1};
    keys{end+1} = fn;
  end;

 case 'clear'
  
  if nargin < 2 % clear all entries
    keys = {};
    data = {};
    
  else % clear entries matching with filename
    i = [];
    for j=1:length(keys)
      if length(keys{j})==length(fn) & ...
	    isequal(keys{j},fn)
	i = j;
      end;
    end;
    
    j = ones(1,length(keys))
    j(i) = 0;
    jj = find(j);
    keys = keys(jj);
    data = data(jj);
    res = [];
  end;
  
 case 'list'
  res = keys;

 case 'limit'
  limsize = fn;
  
 otherwise
  error('command unknown')
end;

s = whos('data');
su = sum(s.bytes);
if su>limsize
  disp(['warning: cache size ',num2str(su),' bytes exceeds limit. Delete entries or rise limit.']);
end;

%| \docupdate 
