classdef onedgrid < gridbase
  % a one dimensional grid implementation

  properties

    global_eind; % global enumeration of entity indices '[1:nelements]'

  end

  methods

    function grid=onedgrid(varargin)
      %function onedgrid(varargin)
      %
      % constructor of a 1d grid
      %
      %     required fields of params:
      %         xnumintervals : number of elements along x directions
      %         xrange : interval covered along the x-axes

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


      % Bernard Haasdonk 18.6.2010

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % copy constructor
      if (nargin>0) & ...
            isa(varargin{1},'onedgrid')
        % copy constructor
        fnames = fieldnames(varargin{1});
        for i=1:length(fnames)
          grid.(fnames{i}) = varargin{1}.(fnames{i});
        end
        % the following only would copy handle!!!
        %grid= varargin{1};
      else
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % default constructor: unit interval
        if (nargin==0)
          params.xrange = [0,1]; % 2 points
          params.xnumintervals = 9;
        else
          params = varargin{1};
        end;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % construct from params

        %  if ~isfield(params,'verbose')
        %    params.verbose = 0;
        %  end;
        %grid = [];
        dx = (params.xrange(2)-params.xrange(1))/params.xnumintervals;
        grid.X = params.xrange(1):dx:params.xrange(2);
        grid.nelements = length(grid.X);
        grid.NBI = [2:grid.nelements, -1; -1,1:(grid.nelements-1)]';
        grid.global_eind = 1:grid.nelements;

      %  grid = class(grid,'onedgrid'); 

      end;
    end

    gridp=gridpart(grid, eind);

    function gcopy = copy(grid);
      % deep copies the grid
      %
      % Return values:
      %   'gcopy': a copy of type onedgrid.
      gcopy=onedgrid(grid);
    end

  end   % methods

end   % classdef
