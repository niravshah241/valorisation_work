classdef ldgdiscfunc < handle
  % an ldg shape functions

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


  properties
    nelements;

    pdeg;

    ndofs_per_element;

    ndofs;

    dimrange;

    dofs;

    grid;
  end

  methods
    function df = ldgdiscfunc(varargin)
      %function df = ldgdiscfunc(varargin)
      % initialize ldg function on triangular grids with input
      %
      % note, the only use of this class is, that by local storage of the
      % dof vector and a (..) as evaluation routine, these objects can be
      % used identically as analytical functions in integration, matrix
      % assembly, etc. But in general the ../ldg directory contains methods
      % for handling ldg functions based on seperate dof and parameter
      % storage. These methods are more efficient as the class&methods.
      %
      % arguments:
      %  varargin: This can be one of
      %   - 'inp = ldgdiscfunc' : copy constructor
      %   - 'inp = {dofs, params}' : 
      %                 optional vector of dofs to be stored in the ldg function
      %   - 'inp = params': 
      %                 Generate an empty ldg function,
      %                 where the required fields of params are:
      %                  - 'pdeg' - polynomial degree on each element
      %                  - 'dimrange' - dimension of function
      %                  - 'nelements' - number of triangles
      %
      % 'df' has a field dofs, which are sorted as follows:
      %  - for all elements
      %    - for all degrees
      %      - for all dimensions
      %
      % i.e. dofs with number '1,1+dimrange,1+2 dimrange' ... are the dofs
      % of the first scalar component, etc.
      %
      % let `\hat \phi_i, i=1...m` be an orthonormal basis on the reference triangle
      % `\hat T`. Let T be an arbitrary triangle and `F_T` be the reference
      % mapping from `\hat T` to `T`. Then for all global dof indices
      % `j=1,...,N` there exists an element `T(j)` and local index `i(j)` such that
      % `\phi_j (x) = \hat \phi_i(j) ( F_T^-1(x))`
      %
      % Then an ldg-discrete function is given by
      %
      % ``
      %  df (x) = \sum_{j=1}^N \text{dof}(j) * \phi_j(x)=
      %           \sum_{j=1}^N \text{dof}(j) * \hat \phi_{i(j)} (F_T(j)^{-1} (x) )
      % ``
      %

      % Bernard Haasdonk 27.1.2009

      if nargin == 0;
        error('no default ldgdiscfunc')
      end;

      if isa(varargin{1},'ldgdiscfunc') % copy constructor

        fnames = fieldnames(varargin{1});
        for i=1:length(fnames)
          df.(fnames{i}) = varargin{1}.(fnames{i});
        end

      else % assume inp argument to be "params" structure

        if nargin == 2
          params = varargin{2};
          dofs = varargin{1};
        else
          params = varargin{1};
          dofs = [];
        end;

        df.nelements = params.nelements;
        df.pdeg = params.pdeg;
        df.ndofs_per_element = ldg_ndofs_per_element(params);
        df.ndofs = ldg_ndofs(params);
        df.dimrange = params.dimrange;

        if isempty(dofs)
          dofs =  zeros(df.ndofs,1);
        end;
        df.dofs = dofs;
        %  keyboard;
      end;

    end

    display(df);

    res = evaluate(df, eindices, lcoord, params);

    res = subsasgn(df, S, val);

    res = subsref(df, S);

  end

end

%| \docupdate
