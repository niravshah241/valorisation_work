classdef cubegrid < gridbase
% a hierarchical cubegrid of arbitrary dimension
%
% This is an axis parallel hexaedral grid, with non-conform refinement
% structure.  All vertices are explicitly generated and stored. Only practical
% for low dimensions, i.e. 1 to 4 or perhaps 5.


  properties
    dimension;        % dimension of grid/number of coordinates per vertex

    vertex;           % 'vertex(i,j)' is the j-th coordinate of i-th vertex point

   % 'vertexindex(id,j)' is the index of j-th vertex of element 'id'
    vertexindex;

    level;            % 'level(id)' equals refinement level of element 'id'

    isleaf;           % 'isleaf(id)' indicates, wether element 'id' is leaf or not

    firstchild;       % 'firstchild(id)' indicates the lowest-index of a child element

    % number of refinement-steps, that have been performed on the grid so far
    refine_steps;

    % 'creation_step(id)' stores the refinement step number that has led to the
    % element 'id'.
    creation_step;

  end

  methods

    function grid= cubegrid(varargin)
    %function cubegrid(varargin)
    % constructor of a hierarchical cubegrid
    %
    % Preliminaries:
    % In the following
    %  - 'gid' denotes global element indices, i.e. also covering non-leaf
    %     elements,
    %  - 'lid' denotes indices for leaf indices.
    %  .
    % Conversion can be performed by 'gids = lid2gid(grid,lid)'
    % Conversion can be performed by 'lids = gid2lid(grid,gid)'
    %
    % The constructor has the following
    % Synopsis:
    %
    %  - 'cubegrid()'      : construction of a default cubegrid (2d unit square)
    %  - 'cubegrid(cgrid)' : copy-constructor
    %  - 'cubegrid(params)': generate cubegrid with certain options. The
    %  argument 'params' requires the following fields:
    %     -# 'range' : cell array of intervals, where the array lengths
    %        determine the grid dimension
    %     -# 'numintervals': vector with number of intervals per dimension
    %
    % perhaps later: constructor by duneDGF-file?
    %
    % generated fields of grid:
    %    dimension     : dimension of grid/number of coordinates per vertex
    %    nelements     : number of overall elements (leaf + nonleaf);
    %    nvertices     : number of vertices;
    %    vertex        : 'vertex(i,j) = 'j-th coordinate of 'i'-th vertex point
    %    vertexindex   : 'vertexindex(id,j) = ' index of 'j'-th vertex of element 'id'
    %    level         : 'level(id) = ' level, i.e. refinement steps of element id
    %    isleaf        : 'isleaf(id) = ' 0 or 1 whether element id is leaf or not
    %    refine_steps  : number of refinement-steps, that have been performed
    %                    on the grid so far
    %    creation_step : 'creation_step(id) ' stores the refinement step number
    %                    that has led to the element 'id'.

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


    % Bernard Haasdonk 1.3.2007

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % default constructor: unit square
      if (nargin==0)
        grid.dimension = 2;
        %range = {[0,1],[0,1]};
        %numintervals = {1,1};
        grid.nelements = 1;
        grid.nvertices = 4;
        % list of vertex vectors: rowwise coordinates
        % => vectex(i,j) : j-th coordinate of i-th vector
        grid.vertex = [0 1 0 1; ...
                       0 0 1 1]';
        % vertex-index vectors: i-th row contains indices of element i
        % element_vertices(i,j): index of j-th vertex of element i
        grid.vertexindex = [1 2 3 4];
        grid.firstchild = 0;
        grid.level = zeros(grid.nelements,1); % everything level 0 by default
        grid.isleaf = ones(grid.nelements,1); % everything leaf by default
        grid.creation_step = zeros(grid.nelements,1); % 0 by default
        grid.refine_steps = 0; % initially set ref-counter to zero

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % copy constructor
      elseif isa(varargin{1},'cubegrid')
        grid= varargin{1};


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % construct from params
      else
        if isnumeric(varargin{1}) && size(varargin{1},1)==1
          partition = varargin{1};
          % creat t partition
          params.numintervals = length(partition)-1;
          params.range = {partition([1,end])};
        else
          partition = [];

          params = varargin{1};
        end
        grid.dimension = length(params.range);
        % abbreviation
        dim = grid.dimension;
        grid.nelements = prod(params.numintervals);
        grid.nvertices = prod(params.numintervals+1);

        % generate vertices
        rangevec ={};
        for i=1:dim
          ni = params.numintervals(i);
          ra = params.range{i};
          rangevec{i} = (0:ni) * (ra(2) - ra(1)) / ni + ra(1);
        end;
        vi = cell(1,dim);
        if (dim > 1)
        [vi{:}] = ndgrid(rangevec{:});
        else
          vi{1} = rangevec{1};
        end;
        ve = zeros(grid.nvertices,0);
        for i=1:dim
          ve = [ve, vi{i}(:)];
        end;
        grid.vertex = ve;

        % generate element vertex index lists
        % list for first element is clear, others obtained by shifting
        % suitably:
        % 1D: [0, 1]
        % 2D: [0, 1, intervals(1)+1, intervals(1)+2]
        % 3D: [ 2D, 2D + (intervals(1)+1)*(intervals(2)+1)],
        % ... doubling plus npoints-product

        ref_elem_indices = [0 1];
        npoints = params.numintervals + 1;
        for i=2:dim
          ref_elem_indices = [ref_elem_indices, ...
                          ref_elem_indices + prod(npoints(1:(i-1)))];
        end;

        % generate list of starting points
        % volume of starting indices:
        % generate vertices
        indexvec ={};
        for i=1:dim
          indexvec{i} = 0:(params.numintervals(i)-1);
        end;
        vi = cell(1,dim);
        if dim > 1
          [vi{:}] = ndgrid(indexvec{:});
        else
          vi{1} = indexvec{1};
        end;
        % starting indices = 1+ vi{1} + npoints(1)*vi{2} + npoints(1)*npoints(2)*vi{3}
        ve = ones(grid.nelements,1);
        for i=1:dim
          ve = ve + prod(npoints(1:(i-1))) * vi{i}(:);
        end;

        % for each starting point: attach shifted array
        grid.vertexindex = ...
        repmat(ve,1,length(ref_elem_indices)) + ...
        repmat(ref_elem_indices,length(ve),1);
        grid.firstchild = zeros(grid.nelements,1);

        grid.level = zeros(grid.nelements,1); % everything level 0 by default
        grid.isleaf = ones(grid.nelements,1); % everything leaf by default
        grid.creation_step = zeros(grid.nelements,1); % 0 by default
        grid.refine_steps = 0; % initially set ref-counter to zero

        if ~isempty(partition)
          grid.vertex(grid.vertexindex(1:end-1,2)) = partition(2:end-1);
        end

      end

    end

    function [compress,new] = tpart_refine(this, lids, new_midpoints)
      assert(this.dimension == 1);
      gids = lid2gid(this, lids);
      refine(this, lids);
      children = this.firstchild(gids);

      if nargin == 3 && ~isempty(new_midpoints)
        this.vertex(this.vertexindex(children,2)) = new_midpoints;
      end
      compress = lids;
      new      = [children,children+1]';
    end


    res = check_consistency(grid);

    leaf_element = coord2leaf_element(grid,coord);

    demo(dummy);

    display(grid);

    ret = get(grid,propertyname, varargin);

    [p0, p1] = get_edges(grid,gid);

    gids = get_leafgids(grid);

    range = get_ranges_of_element(grid, element);

    vols = get_volume(grid,gids);

    lids = gid2lid(grid,gids);

    gids = lid2gid(grid,lids);

    p = plot(grid,params);

    p = plot_grid(grid,params);

    p = plot_leafelement_data(grid,data,params);

    p = plot_leafvertex_data(grid,data,params);

    ngrid = refine(grid, lids);

    ngrid = remove_duplicate_vertices(grid , epsilon );

    function gcopy=copy(grid)
      % a deep copy of the cubegrid
      %
      % return values:
      %   gcopy: a deep copy of this instance also of type cubegrid.
      gcopy = cubegrid(grid);
    end


  end

  methods ( Access = private )

    function gridtmp = gen_plot_data(grid)
      %function gridtmp = gen_plot_data(grid)
      %
      % Auxiliary function for 2d plotting.
      % generation of a structure gridtmp with fields
      % nelements, nvertices, nneigh, VI, X, Y, CX and CY representing
      % the leaf level of the grid
      % as these quantities are required by the common 2d grid plot routines.

      % Bernard Haasdonk 9.5.2007

      elid = get(grid,'leafelements');

      gridtmp = cubegrid();

      gridtmp.nelements = length(elid);
      gridtmp.nvertices = grid.nvertices;
      gridtmp.nneigh = 4;
      gridtmp.VI = grid.vertexindex(elid,[1 2 4 3]);
      gridtmp.X = grid.vertex(:,1);
      gridtmp.Y = grid.vertex(:,2);

      % get X coordinates of vertices as vector
      XX = grid.vertex(gridtmp.VI(:),1);
      % reshape to matrix
      XX = reshape(XX,size(gridtmp.VI));

      % get Y coordinates of vertices as vector
      YY = grid.vertex(gridtmp.VI(:),2);
      % reshape to matrix
      YY = reshape(YY,size(gridtmp.VI));

      gridtmp.CX = mean(XX,2);
      gridtmp.CY = mean(YY,2);
    end


  end
end

