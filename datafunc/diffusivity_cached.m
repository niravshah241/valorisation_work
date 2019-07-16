function [diff] = diffusivity_cached(glob,params,callerid)
%function [diff] = diffusivity_cached(glob,params,callerid)
% compute diffusivity tensor for geometry transformation and store results in a
% cache

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


persistent Dcache hashes cache_filled callerids;
warning('off', 'diffcached:cache:redundance');

%% force clearance of cache
if nargin == 0
  Dcache    = {};
  hashes    = {};
  callerids = {};
  return;
end
% glob column check
if params.debug
  if ~isempty(glob) && size(glob,1) < size(glob,2)
    warning('coordinates in variable glob are given row-wise, but expected them to be column-wise');
    if params.debug > 2
      keyboard;
    end
  end
end

X = glob(:,1);
Y = glob(:,2);

%ishashed = @(x)(min (cellfun(@(y)(all(y==x)), hashes)) == 1);

hasharray = [size(X), X(1), Y(1),callerid];

if params.tstep == 1
  if cache_filled
    Dcache    = {};
    hashes    = {};
    callerids = {};
  end
  cache_filled = false;
  
  [res1, res2] = inv_geo_trans_derivative(params,glob,{(1), (2)},{(1), (2)},callerid);
  row1 = [res1{1}, res1{2}];
  row2 = [res2{1}, res2{2}];

  vlen = size(row1,1);
  temp0  = reshape([ sum(row1 .* row1, 2), sum(row2 .* row2, 2) ]',2*vlen,1);
  tempm1 = reshape([ sum(row2 .* row1, 2), zeros(vlen,1) ]', 2*vlen, 1);
  tempp1 = reshape([ zeros(vlen,1), sum(row1 .* row2, 2) ]', 2*vlen, 1);
  %  temp1 = [ sum(row1 .* row1, 2), sum(row1 .* row2, 2) ];
  %  temp2 = [ sum(row2 .* row1, 2), sum(row2 .* row2, 2) ];
  diff  = spdiags([tempm1,temp0,tempp1],-1:1,2*vlen,2*vlen);

  if(~isempty(hashes))
    hashind = gethash(hasharray, hashes);
  else
    hashind = [];
  end
  if(~(isempty(hashind)) && callerid == callerids{hashind})
    warning('diffcached:cache:redundance', 'two identical hashvalues');
    if(max(max(Dcache{hashind} ~= diff)) ==1 )
        error('WARNING: hashfunction in diffusivity_cached is not good enough!!!');
    end
  else
    % fill cache
    hashind            = length(Dcache)+1;
    hashes{hashind}    = hasharray;
    Dcache{hashind}    = diff;
    callerids{hashind} = callerid;
  end
else
  cache_filled = true;
  hashind = gethash(hasharray, hashes);
  diff = Dcache{hashind};
end

end

function [ind]=gethash(X,hashes)
% function [ind]=gethash(X,hashes)
% compute a hash for the cache function in diffusivity_cached()
ind = find(cellfun(@(y)(all(y==X)), hashes),1);
end

%| \docupdate 
