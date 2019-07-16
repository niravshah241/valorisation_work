function res = quadrature(weights,points,func,varargin)
%function res = quadrature(weights,points,func,varargin)
% 
% integration of function func by given quadrature. func 
% is a function getting a local coordinate vector and varargin
% and giving a (vectorial or scalar) result or a cell array of such.
% points is expected to be a npoints x dimpoint matrix
% result is the same size of f-evaluations.
% we emphasize, that also cell-array valued functions can be
% integrated (for use in rb-methods)

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


% Bernard Haasdonk 27.8.2009

npoints = length(weights);
if isempty(varargin)
  f = func(points(1,:));
  if ~iscell(f)
    res = weights(1)*f;
    for qp=2:npoints
      f = func(points(qp,:));
      res = res +weights(qp)*f;
    end;  
  else % iscell!!
    res = f;
    for q = 1:length(f(:));
      res{q} = weights(1)*f{q};
    end;
    for qp=2:npoints
      f = func(points(qp,:));
      for q = 1:length(f(:));
	res{q} = res{q} +weights(qp)*f{q};
      end;
    end;  
  end;
else % not isempty varargin: same but different...
  f = func(points(1,:),varargin{:});
  if ~iscell(f)
    res = weights(1)*f;
    for qp=2:npoints
      f = func(points(qp,:),varargin{:});
      res = res +weights(qp)*func(points(qp,:),varargin{:});
    end;
  else % iscell
    for q = 1:length(f(:));
      res{q} = weights(1)*f{q};
    end;
    for qp=2:npoints
      f = func(points(qp,:),varargin{:});
      for q = 1:length(f(:));
	res{q} = res{q} +weights(qp)*f{q};
      end;
    end;
  end;
end;

%| \docupdate 
