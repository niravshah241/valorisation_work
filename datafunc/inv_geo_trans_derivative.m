function [P1res, P2res] = inv_geo_trans_derivative(model,glob,P1derivates,P2derivates,callerid)
% function [P1res, P2res] = inv_geo_trans_derivative(model,glob,P1derivates,P2derivates,callerid)
% computes entries of a geometry transformation function's inverse transposed
% jacobian
%
% If 'model.geometry_transformation' is set to 'spline', this function computes
% derivatives of the given geometry transformation function 
% `\Phi : \mathbb{R}^2 \to \mathbb{R}^2`.
%
% Parameters:
%  glob:       matrix with two columns for the `x` and `y` coordinates
%              specifying where the derivative should be evaluated.
%  P1derivates: a cell array of vectors of length one to three filled with
%              scalar values one or two. For each cell entry, a derivative is
%              evaluated, where a value of one means a derivation in `x`
%              direction and two a derivation in `y` direction. For example
%              '{ [1, 2, 1] }' corresponds to the derivative
%              `\partial_x \partial_y \partial_x \Phi_1(x,y)`
%  P2derivates: a cell array of vectors of length one to three filled with
%              scalar values one or two. For each cell entry, a derivative is
%              evaluated, where a value of one means a derivation in `x`
%              direction and two a derivation in `y` direction. For example
%              '{ [1, 2, 1] }' corresponds to the derivative
%              `\partial_x \partial_y \partial_x \Phi_2(x,y)`
%  callerid:   As this function only depends on the space variable, during the
%              evaluation of an evolution scheme, the same derivatives need to
%              be computed repeatedly. Therefore, the computations done in the
%              first time step are are cached by this function. In order to
%              make sure, the hash function works correctly, all calls of this
%              function with different arguments for 'P1derivates' and
%              'P2derivates' should pass a unique 'callerid'.
%
% Return values:
%  P1res:      a cell array of derivative evaluations specified by 'glob' and
%              'P1derivates'.
%  P2res:      a cell array of derivative evaluations specified by 'glob' and
%              'P2derivates'.

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


persistent P1cache P2cache hashes cache_filled callerids;
warning('off','geotrans:cache:redundance');

%% force clearance of cache
if nargin == 0
  P1cache = {};
  P2cache   = {};
  hashes    = {};
  callerids = {};
  return;
end

%ishashed = @(x)(min (cellfun(@(y)(all(y==x)), hashes)) == 1);
%gethash  = @(x)(find(cellfun(@(y)(all(y==x)), hashes)));

%cleanstring = @(x)(strrep(strrep(strrep(strrep(x,'-','m'),'.','p'),' ',''),'+','P'));

% ignore ID and just compute a hash function (size and first coordinate pair should be enough)
%hashvalue = ['h', strrep(strrep(sprintf('%dx%d', size(X)),'.','p'),'-','m'), ...
%                  strrep(strrep(sprintf('a%dx%d',X(1),Y(1)),'.','p'),'-','m')];
X = glob(:,1);
Y = glob(:,2);
hasharray = [size(X), X(1), Y(1), P1derivates{1:2}, P2derivates{1:2}, length(P1derivates),callerid];
%hashvalue = ['h', cleanstring(sprintf('%dx%d', size(X))), ...
%                  cleanstring(sprintf('a%dx%d',X(1),Y(1))), ...
%                  cleanstring(num2str([P1derivates{:}])), ...
%                  cleanstring(num2str([P2derivates{:}]))];

if model.tstep == 1
%  if isfield(P1cache, 'cache_is_filled') && P1cache.cache_is_filled
%    P1cache = rmfield(P1cache,fieldnames(P1cache));
%    P2cache = rmfield(P2cache,fieldnames(P2cache));
%  end
  if cache_filled
    if model.verbose > 19
      disp('erase cache values');
    end
    P1cache   = {};
    P2cache   = {};
    hashes    = {};
    callerids = {};
  end
  cache_filled = false;
  if isequal(model.geometry_transformation, 'spline')

    %    x_hill  = model.geometry_transformation_spline_x;
    %    y_hill  = model.geometry_transformation_spline_y;
    %    y_hill(2) = model.hill_height;
    p_mu    = spline_select(model);
    [ breaks, coeffs, pieces, order ] = unmkpp(p_mu);
    p_mu_d  = mkpp(breaks, coeffs(:,1:order-1) .* repmat(order-1:-1:1,pieces,1));
    [ breaks, coeffs, pieces, order ] = unmkpp(p_mu_d);
    if order == 1
      p_mu_dd = mkpp(breaks, zeros(size(coeffs)));
    else
      p_mu_dd = mkpp(breaks, coeffs(:,1:order-1) .* repmat(order-1:-1:1,pieces,1));
    end

    phi_d = cell(2,2);
    phi_d{1,1} = @(X,Y) (ones(size(X)));
    phi_d{1,2} = @(X,Y) (zeros(size(X)));
    phi_d{2,1} = @(X,Y) (-Y.*ppval(p_mu_d, X)./(1 + ppval(p_mu, X)));
    phi_d{2,2} = @(X,Y) (1./(1+ppval(p_mu, X)));

    phi_dd = cell(2,2,2);
    phi_dd{1,1,1} = @(X,Y) (zeros(size(X)));
    phi_dd{1,1,2} = @(X,Y) (zeros(size(X)));
    phi_dd{2,2,1} = @(X,Y) -ppval(p_mu_d, X)./(1 + ppval(p_mu, X));
    phi_dd{2,2,2} = @(X,Y) (zeros(size(X)));
    phi_dd{1,2,1} = @(X,Y) (zeros(size(X)));
    phi_dd{2,1,2} = @(X,Y) -ppval(p_mu_d, X)./(1 + ppval(p_mu, X));

    phi_ddd = cell(2,2,2,2);
    phi_ddd{1,1,1,1} = @(X,Y) (zeros(size(X)));
    phi_ddd{1,1,2,1} = @(X,Y) (zeros(size(X)));
    phi_ddd{1,2,1,1} = @(X,Y) (zeros(size(X)));
    phi_ddd{1,2,2,1} = @(X,Y) (zeros(size(X)));
    phi_ddd{2,1,1,2} = @(X,Y) (ppval(p_mu_dd, X)./(1+ppval(p_mu, X)) - ppval(p_mu_dd, X) ./ (1+ppval(p_mu, X) ));
    phi_ddd{2,1,2,2} = @(X,Y) (zeros(size(X)));
    phi_ddd{2,2,1,2} = @(X,Y) (zeros(size(X)));
    phi_ddd{2,2,2,2} = @(X,Y) (zeros(size(X)));

    P1res = cell(size(P1derivates));
    for i = 1 : length(P1derivates)
      p1i     = P1derivates{i};
      lenderi = length(p1i);
      if lenderi == 1
        P1res{i} = phi_d{1,p1i(1)}(X,Y);
      elseif lenderi == 2
        P1res{i} = phi_dd{1,p1i(1),p1i(2)}(X,Y);
      else
        P1res{i} = phi_ddd{1,p1i(1),p1i(2), p1i(3)}(X,Y);
      end
    end
    P2res = cell(size(P2derivates));
    for i = 1 : length(P2derivates)
      p2i     = P2derivates{i};
      lenderi = length(p2i);
      if lenderi == 1
        P2res{i} = phi_d{2,p2i}(X,Y);
      elseif lenderi == 2
        P2res{i} = phi_dd{2,p2i(1),p2i(2)}(X,Y);
      else
        P2res{i} = phi_ddd{2,p2i(1),p2i(2), p2i(3)}(X,Y);
      end
    end
    if(~isempty(hashes))
      hashind = gethash(hasharray, hashes);
    else
      hashind = [];
    end
    if(~( isempty(hashind)) && callerid == callerids{hashind})
      warning('geotrans:cache:redundance','two identical hashvalues');
      if(max(max([P1cache{hashind}{:}] ~= [P1res{:}])) == 1 || ...
           max(max([P2cache{hashind}{:}] ~= [P2res{:}])) == 1)
        error('WARNING: hashfunction is not good enough!!!');
      end
    else
      %      if(isfield(P1cache, hashvalue))
      %% fill cache
      hashind = length(P1cache)+1;
      hashes{hashind} =hasharray;
      P1cache{hashind}=P1res;
      P2cache{hashind}=P2res;
      callerids{hashind}=callerid;
    end
  end
else  %  model.tstep != 1
  cache_filled = true;
  hashind = gethash(hasharray, hashes);
  P1res = P1cache{hashind};
  P2res = P2cache{hashind};
  %  P1res = P1cache.(hashvalue);
  %  P2res = P2cache.(hashvalue);
end

end

function [ind]=gethash(Xx,hashes)
% function [ind]=gethash(Xx,hashes)
% private function computing a hash index from the function arguments to
% inv_geo_trans_derivative().
  ind = find(cellfun(@(y)(all(y==Xx)), hashes),1);
end

