function ENBI = get_enbi(grid, edge, tstep)
%function ENBI = get_enbi(grid, edge)
%
% function assembling a matrix with the 5 neighbour's cell indices that
% are needed in order to compute the gradient over the edge given by
% 'edge' in each row. see also the sketch below.
%
% arrangement of cell indices:
%  edges' indices of the main cell by which one gets the "right"
%  neighbours 
% @verbatim
% Example: edge == 1:
%   N   |   NE                        ENBI(:,2) = N  | ENBI(:,5) = NE
% ------|-------   <-- nb_ind(0)     ----------------|-----------------
%  main |   E                         ENBI(:,1)      | ENBI(:,4)
% ------|-------   <-- nb_ind(1)     ----------------|-----------------
%   S   |   SE                        ENBI(:,3)      | ENBI(:,6)
%       ^
%      edge
% @endverbatim
%
% Martin Drohmann 05.03.2008

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


persistent Ecache hashes cache_filled;
warning('off','enbi:cache:redundance');

if nargin == 1
  Ecache = {};
  hashes = {};
  return;
end

%ishashed = @(x)(min (cellfun(@(y)(all(y==x)), hashes)) == 1);
hasharray = [size(grid.X), grid.X(1), grid.Y(1), edge];

if tstep == 1 
  if cache_filled
%    disp('erase cache values');
    Ecache = {};
    hashes  = {};
  end
  cache_filled = false;
  num_cells = size(grid.CX,1);

  % TODO: this is a stupid hack: I don't like it
  temp = edge;
  if edge > 2;
    temp = edge - 2;
  end
  nb_ind = mod([ 2 0 ] + temp, 4) + 1;


  ENBI = zeros(num_cells, 6);
  ENBI(:,1)     = 1:num_cells;
  ENBI(:,[2,3]) = grid.NBI(1:num_cells, nb_ind);

  %neg_ind       = find(ENBI < -1);
  %[ind,col]     = ind2sub(size(ENBI), neg_ind);
  %ENBI(neg_ind) = ind;

  for i = 1:3
    non_bnd_ind = ENBI(:,i) > 0;
    dir_bnd_ind = ENBI(:,i) == -1;
    ENBI(non_bnd_ind,i+3) = grid.NBI(ENBI(non_bnd_ind,i),edge);
    % TODO: Also assign to -1 if _one_ neighbour already is -1. It should be
    % enough to check for index 4
    ENBI(dir_bnd_ind,i+3) = -1;
  end

  temp                    = logical(ENBI(:,4) == -1);
  corner_bnd_ind          = temp & ENBI(:,5) == 0;
  ENBI(corner_bnd_ind, 5) = -1; 
  corner_bnd_ind          = temp & ENBI(:,6) == 0;
  ENBI(corner_bnd_ind, 6) = -1; 
  %neg_ind       = find(ENBI < -1);
  %[ind,col]     = ind2sub(size(ENBI), neg_ind);
  %ENBI(neg_ind) = ENBI(sub2ind(size(ENBI), ind, col - 3));

  %ENBI(find(ENBI == 0)) = -1;

  if(~isempty(hashes))
    hashind = gethash(hasharray, hashes);
  else
    hashind = [];
  end
  if(~( isempty(hashind)))
    warning('enbi:cache:redundance','two identical hashvalues');
    %       if(max(max([P1cache.(hashvalue){:}] ~= [P1res{:}]))==1 || ...
    %           max(max([P2cache.(hashvalue){:}] ~= [P2res{:}]))==1)
    if(max(max(Ecache{hashind} ~= ENBI)) == 1)
      error('WARNING: hashfunction is not good enough!!!');
    end
  else
    hashind = length(Ecache)+1;
    hashes{hashind} = hasharray;
    Ecache{hashind} = ENBI;
  end

else

  cache_filled = true;
  hashind = gethash(hasharray, hashes);
  ENBI = Ecache{hashind};
end
end

function [ind]=gethash(X,hashes)
  ind = find(cellfun(@(y)(all(y==X)), hashes));
end

