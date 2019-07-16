classdef gridbase < handle
  % Base class for all grid classes
  %

  properties

    nelements; % number of overall elements (leaf + nonleaf)

    nvertices; % number of vertices

    nneigh; % number of neighbor elements

    A; % vector of element area

    Ainv; % vector of inverted element area

    % matrix of vertex indices: 'VI(i,j)' is the global index of
    % 'j'-th vertex of element 'i'
    VI;

    X; % vector of vertex x-coordinates

    Y; % vector of vertex y-coordinates

    CX; % vector of centroid x-values

    CY; % vector of centroid y-values

    % 'NBI(i,j)' = element index of 'j'-th neighbour of element 'i'
    % boundary faces are set to -1 or negative values are requested
    % by 'params.boundary_type'
    NBI;

    % 'INB(i,j)' = local edge number in 'NBI(i,j)' leading from
    % element 'NBI(i,j)' to element 'i', i.e. satisfying
    % 'NBI(NBI(i,j), INB(i,j)) = i'
    INB;

    % 'EL(i,j)' = length of edge from element 'i' to neighbour 'j'
    EL;

    % 'DC(i,j)' = distance from centroid of element 'i' to NB 'j'
    % for boundary elements, this is the distance to the reflected
    % element (for use in boundary treatment)
    DC;

    % 'NX(i,j)' = x-coordinate of unit outer normal of edge from el
    % 'i' to NB 'j'
    NX;

    % 'NY(i,j)' = y-coordinate of unit outer normal of edge from el
    % 'i' to NB 'j'
    NY;

    % 'ECX(i,j)' = x-coordinate of midpoint of edge from el 'i' to
    % NB 'j'
    ECX;

    % 'ECY(i,j)' = x-coordinate of midpoint of edge from el 'i' to
    % NB 'j'
    ECY;

    % vector of x-coordinates of point `S_i` (for rect: identical to
    % centroids)
    %
    % for diffusion-discretization with FV-schemes, points `S_i`
    % must exist, such that `S_i S_j` is perpendicular to edge `i j`
    % the intersection are denoted `S_ij`
    SX;

    % vector of y-coordinates of point `S_i` (for rect: identical to
    % centroids)
    %
    % for diffusion-discretization with FV-schemes, points `S_i`
    % must exist, such that `S_i S_j` is perpendicular to edge `i j`
    % the intersection are denoted `S_ij`
    SY;

    % 'ESX(i,j)' = x-coordinate of point `S_ij` on edge el 'i' to
    % NB 'j'
    ESX;

    % 'ESY(i,j)' = x-coordinate of point `S_ij` on edge el 'i' to
    % NB 'j'
    ESY;

    % 'DS(i,j)' = distance from `S_i` of element 'i' to `S_j` of NB
    % 'j' for boundary elements, this is the distance to the
    % reflected element (for use in boundary treatment)
    DS;

    hmin; % minimal element-diameter

    % geometry bound (simultaneously satisfying
    % ``\alpha h_i^d \leq A(T_i),``
    % ``\alpha \text{diam}(T_i) \leq h_i^{d-1}`` and
    % ``\alpha h_i \leq \text{distance(midpoint i to any neighbour)}``
    alpha;

    % Jacobian inverse transposed 'JIT(i,:,:)' is a '2x2'-matrix of
    % the Jacobian Inverse Transposed on element 'i'
    JIT;

  end

  methods

    function obj = gridbase
      % default constructor for gridbase

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

    end


    res = check_consistency(grid);

    display(grid);

    F = edge_quad_eval(grid,elids,edgeids,degree,FF);

    F = edge_quad_eval_mean(grid,elids,edgeids,degree,FF);

    PP = edge_quad_points(grid,elids,edgeids,degree);

    [P1, P2] = get_edge_points(grid,elids,edgeids);

    ENBI = get_enbi(grid, edge, tstep);

    gridp=gridpart(grid, eind);

    inspect(grid, params);

    p = plot_polygon_grid(grid,params);

    p = plot_element_data(grid,data,plot_params);

    plot_element_data_sequence(grid,data,plot_params);

    p = plot_vertex_data(grid,data,params);

    grid = set_boundary_types(grid,params);

  end

  methods(Abstract)

    gcopy = copy(grid);
    % returns a deep copy object of the grid implementation
    %
    % Every grid implementations needs to provide a deep copy method, because
    % it is derived from a handle class.
    %
    % Return values:
    %   gcopy: Deep copy of type gridbase.

  end

end
