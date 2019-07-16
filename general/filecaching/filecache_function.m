function [varargout] = filecache_function(funcptr, varargin)
%function [varargout] = filecache_function(funcptr, varargin)
% function used for file-caching other function calls.
%
% If an expensive function
%
% @code
%     [E, F] = myfunction(A,B,C)
% @endcode
%
% is called frequently during a program run, a filecaching can be
% used, i.e. one calls the function instead as
%
% @code
%     [E, F] = filecache_function(@myfunction,A,B,C);
% @endcode
%
% If the function result exists in the cache, this is loaded, otherwise the
% function is called, the result saved in the cache and the result returned.
%
% Parameters:
%  funcptr:  is the pointer to a function whose calls are to be cached.
%  varargin: is parameter list of the cached function call
%
% Return values:
%  varargout: return values of the cached function call
%
% - The cache-directory is assumed to be in 'RBMATLABTEMP/cache'
% - The cache-directory can be cleared with the function filecache_clear()

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


% Bernard Haasdonk 22.5.2007

if nargout==0
  error('filecaching only works for functions with return values!');
end;

funcname = func2str(funcptr);

% determining of a 'hash-code' from the function-arguments:
tmpfn = fullfile(filecache_path,'tmparg.mat');
saveargs = varargin;
for i=1:length(saveargs)
  if isobject(varargin{i})
    warning('off', 'MATLAB:structOnObject');
    saveargs{i} = struct(varargin{i});
    warning('on', 'MATLAB:structOnObject');
  end
end

save(tmpfn,'saveargs');
fid = fopen(tmpfn,'r');
bin = fread(fid);
% cut of header, which contains date, i.e. first 116 bytes
bin = bin(117:end);
fclose(fid);
uintbin = uint32(bin);
% simple arithmetics for identifying arguments 
key  = mod(sum(uintbin.*uint32(1:length(uintbin))'),1e9);
keystr = num2str(key);
%keyboard;
%key = sum()

resultfn      = fullfile(filecache_path,...
                         [funcname,keystr,'.mat']);
multval_cache = fullfile(filecache_path,...
                         [funcname,'mv.mat']);

is_cached  = 0;
hash_found = 0;
funcstr = strtrim(evalc('disp(funcptr)'));

if exist(resultfn,'file');
  loaded = load(resultfn);
  hash_found = 1;
elseif exist(multval_cache,'file')
  tmp = load(multval_cache);
  tmp_key_pos = find(tmp.keylist==key,1);
  if ~isempty(tmp_key_pos)
    tmp_key_ind  = tmp.map_key2ind(tmp_key_pos);
    tmp_file_ind = tmp.map_key2file(tmp_key_pos);
    storefile = fullfile(filecache_path,...
                         [funcname,'mv',num2str(tmp_file_ind),'.mat']);
    if exist(storefile,'file')
      tmp = load(storefile);
      loaded = tmp.cache(tmp_key_ind);
      hash_found = 1;
    else
      error(['Could not find file ',storefile,...
            ', which is specified in multivalue cache list.']);
    end
  end
end

if hash_found
  % isequal does not work for handles
  %  if isequal(tmp.saveargs,saveargs)
  if isequal_ignore_handles(loaded.saveargs,saveargs)
    varargout = loaded.varargout;
    is_cached = 1;
    disp(['call of ',funcstr,', successfully read from cache, key=',keystr]);
%    keyboard;
  else
    % the key-computation formula must be updated in this case...
    disp(['arguments in cached file and current call do not',...
          'correspond!!']);
    disp('please change key-function in file_cache_function!');
    disp('recomputation is started and cache-file replaced.');
    %keyboard;
  end
end


if ~is_cached
  % call function
  disp('result not found in cache');
  [varargout{1:nargout}] = funcptr(varargin{:});
  %  save(resultfn,'varargin');

  siz=whos('varargout');
  % if size is less than 100kb, store in multi-value cache.
  if(siz.bytes < 1024*100)
    if exist(multval_cache,'file')
      meta=load(multval_cache);
      max_file_ind = meta.map_key2file(end);
      max_key_ind  = meta.map_key2ind(end)+1;
      last_size    = meta.last_size;
      % if store file size exceeds 10MB, create a new one
      if (last_size + siz.bytes > 1024*10000)
        max_file_ind = max_file_ind + 1;
        max_key_ind  = 1;
      end
    else
      meta.keylist      = [];
      meta.map_key2file = [];
      meta.map_key2ind  = [];
      max_file_ind      = 1;
      max_key_ind       = 1;
    end
    storefile = fullfile(filecache_path,...
                         [funcname,'mv',num2str(max_file_ind),'.mat']);
    if exist(storefile, 'file')
      load(storefile);
    end
    cache(max_key_ind).saveargs  = saveargs;
    cache(max_key_ind).varargout = varargout;

    tmp_size  = whos('cache');
    last_size    = tmp_size.bytes;
    keylist      = [meta.keylist, key];
    map_key2file = [meta.map_key2file, max_file_ind];
    map_key2ind  = [meta.map_key2ind, max_key_ind];

    save(multval_cache, 'last_size', 'keylist', 'map_key2file', 'map_key2ind');
    save(storefile, 'cache');

  else
    save(resultfn,'saveargs','varargout');
  end
  disp(['call of ',funcstr,', now computed and cached, key=',keystr]);
end

