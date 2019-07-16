classdef femdiscfunc < handle
% class femdiscfunc
% class representing a continous piecewise polynomial function of arbitrary
% dimension. 'DOFS' correspond to the values of Lagrange-nodes.
% 
% \link femdiscfunc::global_dof_index 'global_dof_index(elid,lagrange_node)'
% \endlink yields the global index of the first dof of the basis function
% corresponding to the given Lagrange node and element elid.
% 'gid:(gid+dimrange-1)' are the subsequent dofs of the vectorial function in
% the lagrange node.  first all dofs in nodes are counted, then all dofs in
% element interior, then the dofs on edge-interiors.
%
% The Lagrange nodes 'l_1,...,l_m' with 'm=0.5*(pdeg+1)*(pdeg+2)'
% are sorted in the following order:
% @verbatim
%       l_m = v_3
%       *
%       |\
%       | \
%       |  \
%       *   *
%       |    \
%       |     \
%       |______\
%       *   *   *
% v_1 = l_1      v_2 = l_(pdeg+1)
% @endverbatim
%
% where 'v_1, v_2, v_3' denote the sorting of the triangles corners.

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


% Bernard Haasdonk 11.1.2011

  properties
    pdeg;
    dimrange;
    grid;
    df_info;
    dofs;
    % dependent variables:

    % nelements documentation
    nelements;
    ndofs; 
    ndofs_per_element;
    % the global dof index
    global_dof_index;

  end

  methods

    function df = femdiscfunc(dofs,df_info)
    % constructor, dofs possibly [], grid required!
    %
    % parameters:
    %   dofs:     a vector of #'df.ndofs' global degress of freedom for the
    %             discrete function, if '==[]' a zero vector is created.
    %   df_info:  object of type ::feminfo describing the structure of the
    %             underlying discrete function space
    %
    % required fields of df_info:
    %   pdeg:     polynomial degree of lagrange functions
    %   dimrange: dimension of the range
    %   grid.nelements: @copybrief gridbase::nelements
      df.pdeg = df_info.pdeg;
      df.dimrange = df_info.dimrange;
      df.grid = df_info.grid;
      df.df_info = df_info;
      df.dofs=dofs;
      % dependent:
      df.nelements = df_info.grid.nelements;
      df.ndofs = df_info.ndofs;
      df.ndofs_per_element = df_info.ndofs_per_element;
      if isempty(dofs)
        df.dofs = zeros(df.ndofs,1);
      end;
      df.global_dof_index = df_info.global_dof_index;
    end

    function sdf = scalar_component(df,ncomp)
    % extraction of component of vectorial fem function
      [scalardofs, scalar_df_info] = ...
        fem_scalar_component(df.dofs, ncomp, df);
      sdf = femdiscfunc(scalardofs,scalar_df_info);
    end

    function p = plot(df,params)
    % plot as colormap
      if nargin<2
        params = [];
      end;
    p = plot_discfunc(df,params);
    end;

    function p = plot_dofmap(df,params)
    % plot as colormap
      if nargin<2
        params = [];
      end;
      p = fem_plot_dofmap(df,params);
    end

    function res = evaluate(df,einds,lcoord,dummy1,dummy2)
    % plot as colormap
    %
    % description of evaluate function
      res = fem_evaluate(df,einds,lcoord,[],[]);
    end

    function res = subsref(df, S)
    % This method enables indexation of discrete functions
    %
    % redirects arguments to evalute function, so
    % @code
    %   df(einds, lcoord) == evaluate(df, einds, lcoord)
    % @endcode

      %    t = S.type;
      %    keyboard;
      if S(1).type~='('
        res = builtin('subsref', df, S);
      else
        res = evaluate(df,S.subs{:});
      end;
    end
    
    % copy
    function cdf = copy(df)
    cdf = femdiscfunc(df.dofs,df.df_info);    
    end;
    
    % addition
    function df3 = plus(df1,df2)
    df3 = femdiscfunc([],df1.df_info);
    df3.dofs = df1.dofs + df2.dofs;
    end;

    % subtraction
    function df3 = minus(df1,df2)
    df3 = femdiscfunc([],df1.df_info);
    df3.dofs = df1.dofs - df2.dofs;
    end;

    % negative
    function df2 =uminus(df1)
    df2 = femdiscfunc([],df1.df_info);
    df2.dofs = -df1.dofs;
    end;

    % multiplication
    function df2 =mtimes(factor,df1)
    df2 = femdiscfunc([],df1.df_info);
    df2.dofs = factor * df1.dofs;
    end;

  end  % methods

end % classdef
