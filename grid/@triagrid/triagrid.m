classdef triagrid < gridbase
  % a triangular conforming grid in two dimesions
  %
  % General geometrical information is stored including neighbour
  % information. Therefore also boundary neighbour relations can be
  % specified. Boundary edges can be marked by "painting rectangles": the
  % edges with midpoints within such a rectangle are marked accordingly.
  % By this boundary edges can be marked for later special use.  Much
  % (partially redundant) information is stored in the grid, which might be
  % useful in simulations.

  properties
     nedges_interior;
     nedges_boundary;
  end;

  methods

    function grid= triagrid(varargin)
      %function triagrid(varargin)
      %
      % constructor of a triangular conform grid in 2 dimensions with following
      % Synopsis:
      %
      % - triagrid() : construction of a default triagrid (2d unit square,
      %              loaded from file, -1 as outer neighbour indices)
      % - @link triagrid() triagrid(tgrid) @endlink : copy-constructor
      % - @link triagrid() triagrid(params) @endlink : in this case either
      %    -# the field grid_initfile is existent in params. Then the file is read
      %       and a pointlist p and a triangle. Procedure is then executing the
      %       following constructor
      %    -# a structured equidistant triangle grid is generated
      %        - 'params.xnumintervals' : number of elements along x directions
      %        - 'params.ynumintervals' : number of elements along y directions
      %        - 'params.xrange,yrange' : intervals covered along the axes
      %        .
      %        with the optional fields
      %        - 'params.bnd_rect_corner1' : coordinates of lower corner of to
      %        be marked boundaries
      %        - 'params.bnd_rect_corner2' : coordinates of upper corner of to
      %        be marked boundaries
      %        - 'params.bnd_rect_index': integer index to be set on the edges
      %        in the above  defined rectangle. Should not be positive integer
      %        in the range of the number of elements. use negative indices for
      %        certain later discrimination.
      %        .
      %        For the last three optional boundary settings, also multiple
      %        rectangles can be straightforwardly defined by accepting matrix
      %        of columnwise corners1, corners2 and a vector of indices for the
      %        different rectangles.
      % - @link triagrid() triagrid(p,t,params) @endlink : generate triagrid from
      %   triangle-data with certain options.
      %       - p is assumed to be a 2 x npoints matrix with coordinates
      %       - t is assumed to be a XX x ntriangles matrix, but only first three
      %         rows are used == vertex indices, in clockwise order as default, all
      %         nondefined edges are set to -1, then the following "rectangles" are
      %         set as specified in params
      %
      % Using this class, grids from PDETOOLS can be used:
      %
      %     pdetools
      %     %%% => generate your grid and export 'p' and 't' to MATLAB workspace
      % @code
      %     save('mygrid','p','t')
      %     grid = triagrid(struct('grid_initfile','mygrid'));
      % @endcode
      %
      % perhaps later: constructor by duneDGF-file?
      % perhaps later: contructor-flag: full vs non-full
      %                 => only compute redundant information if required.
      %
      % Parameters:
      %  'varargin' : variable number of input arguments, see above for description
      %  of possible configurations.
      %
      % Return values:
      %  'grid' : generated triagrid
      %
      % generated fields of grid:
      % nelements: number of elements
      % nvertices: number of vertices
      % nneigh: 3
      % A    : vector of element area
      % Ainv : vector of inverted element area
      % X    : vector of vertex x-coordinates
      % Y    : vector of vertex y-coordinates
      % VI   : matrix of vertex indices: 'VI(i,j)' is the global index of j-th
      %        vertex of element i
      % CX   : vector of centroid x-values
      % CY   : vector of centroid y-values
      % NBI  : 'NBI(i,j) = ' element index of j-th neighbour of element i
      %        boundary faces are set to -1 or negative values are requested by
      %        'params.boundary_type'
      % INB  : 'INB(i,j) = ' local edge number in 'NBI(i,j)' leading from
      %        element 'NBI(i,j)' to element 'i', i.e. satisfying
      %        'NBI(NBI(i,j), INB(i,j)) = i'
      % EL   : 'EL(i,j) = ' length of edge from element 'i' to neighbour 'j'
      % DC   : 'DC(i,j) = ' distance from centroid of element i to NB j
      %        for boundary elements, this is the distance to the reflected
      %        element (for use in boundary treatment)
      % NX   : 'NX(i,j) = ' x-coordinate of unit outer normal of edge from el
      %        'i' to NB 'j'
      % NY   : 'NY(i,j) = ' y-coordinate of unit outer normal of edge from el
      %        'i' to NB 'j'
      % ECX  : 'ECX(i,j) = ' x-coordinate of midpoint of edge from el 'i' to NB 'j'
      % ECY  : 'ECY(i,j) = ' y-coordinate of midpoint of edge from el 'i' to NB 'j'
      % SX   : vector of x-coordinates of point `S_i` (for rect: identical to
      %        centroids)
      % SY   : vector of y-coordinate of point `S_j` (for rect: identical to
      %        centroids)
      % ESX  : 'ESX(i,j) = ' x-coordinate of point `S_{ij}` on edge el i to NB j
      % ESY  : 'ESY(i,j) = ' y-coordinate of point `S_{ij}` on edge el i to NB j
      % DS   : 'DS(i,j) = ' distance from `S_i` of element i to `S_j` of NB j
      %        for boundary elements, this is the distance to the reflected
      %        element (for use in boundary treatment)
      % hmin : minimal element-diameter
      % alpha: geometry bound (simultaneously satisfying `\alpha \cdot h_i^d
      %        \leq A(T_i)`, `\alpha \cdot \mbox{diam}(T_i) \leq h_i^{d-1}` and
      %        `\alpha \cdot  h_i \leq `dist(midpoint `i` to any neigbour) )
      % JIT  : Jacobian inverse transposed 'JIT(i,:,:)' is a 2x2-matrix of the Jac.
      %        Inv. Tr. on element 'i'
      %
      % \note for diffusion-discretization with FV-schemes, points `S_i` must exist,
      % such that `S_i S_j` is perpendicular to edge 'i j', the intersection points
      % are denoted `S_{ij}`
      %

% This program is open source.  For license terms, see the COPYING file.
%
% --------------------------------------------------------------------
% ATTRIBUTION NOTICE:
% This product includes software developed for the RBmatlab project at
% (C) Universities of Stuttgart and MÃ¼nster, Germany.
%
% RBmatlab is a MATLAB software package for model reduction with an
% emphasis on Reduced Basis Methods. The project is maintained by
% M. Dihlmann, M. Drohmann, B. Haasdonk, M. Ohlberger and M. Schaefer.
% For Online Documentation and Download we refer to www.morepas.org.
% --------------------------------------------------------------------


      % Bernard Haasdonk 9.5.2007

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % copy constructor
      if (nargin==1) & isa(varargin{1},'triagrid')
        grid = triagrid; % create default triagrid

        fnames = fieldnames(varargin{1});
        for i=1:length(fnames)
          grid.(fnames{i}) = varargin{1}.(fnames{i});
        end

      else

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % default constructor: unit square
        if (nargin==0)
          %    p = [0.0   0.5  1.0   0.0    0.5   1.0   0.0    0.5    1.0; ...
          %	 0.0   0.0  0.0   0.5    0.5   0.5   1.0    1.0    1.0 ];
          %     t = [1   2  2  3  4  5  5  6; ...
          %	  2   5  3  6  5  8  6  9; ...
          %	  4   4  5  5  7  7  8  8] ;
          params = [];
          load 2dsquaretriang;

          % mark boundary in rectangle from [-1,-1] to [+2,+2] with index -1
          params.bnd_rect_corner1 = [-1,-1]';
          params.bnd_rect_corner2 = [2 2]';
          params.bnd_rect_index = [-1];
        elseif nargin==1
          params = varargin{1};
          if isfield(params,'grid_initfile')
            tmp = load(params.grid_initfile);
            t = tmp.t(1:3,:);
            p = tmp.p;
          else
            nx = params.xnumintervals;
            ny = params.ynumintervals;
            nvertices =(nx+1)*(ny+1);

            dx = (params.xrange(2)-params.xrange(1))/nx;
            dy = (params.yrange(2)-params.yrange(1))/ny;

            % set vertex coordinates
            vind = (1:(nx+1)*(ny+1));
            X = mod(vind-1,nx+1) * dx + params.xrange(1);
            Y = floor((vind-1)/(nx+1))*dy  + params.yrange(1);
            p = [X(:)'; Y(:)'];

            % identify triangles:
            % ----
            % | /|
            % |/ |
            % ----

            all_ind = 1:(nx+1)*(ny+1);
            % index vector of lower left corners, i.e. not last col/row
            lc_ind = find((mod(all_ind,nx+1)~=0) & (all_ind < (nx+1)*ny));
            t1 = [lc_ind;lc_ind+1;lc_ind+2+nx]; % lower triangles
            t2 = [lc_ind;lc_ind+nx+2;lc_ind+nx+1]; % upper triangles
            t = [t1,t2];
          end;
        else % read p-e-t information from varargin
          p =varargin{1};
          t =varargin{2};
          if size(t,1)>3
            t = t(1:3,:);
          end;
          params = varargin{3};
        end;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % construct from p,t and params

        if ~isfield(params,'verbose')
          params.verbose = 0;
        end;

      %  nx = params.xnumintervals;
      %  ny = params.ynumintervals;

        grid.nelements = size(t,2);
        grid.nvertices = size(p,2);
        grid.nneigh= 3;

        % set vertex coordinates
        grid.X = p(1,:)';
        grid.Y = p(2,:)';

        % set element vertex indices: numbering counterclockwise
        % counterclockwise, edge j connects vertex j and j+1

        grid.VI = t';

        % get areas of grid-cells after VI and X,Y are available
        %
        %  A(a b c) = 0.5 |Det ((b-a) (c-a)) |
        grid.A = abs(area_triangle(...
            [grid.X(grid.VI(:,1)), ...
             grid.Y(grid.VI(:,1))], ...
            [grid.X(grid.VI(:,2)), ...
             grid.Y(grid.VI(:,2))], ...
            [grid.X(grid.VI(:,3)), ...
             grid.Y(grid.VI(:,3))]));
        grid.A = grid.A(:);
      %  a11 = grid.X(grid.VI(:,2))-grid.X(grid.VI(:,1));
      %  a21 = grid.Y(grid.VI(:,2))-grid.Y(grid.VI(:,1));
      %  a12 = grid.X(grid.VI(:,3))-grid.X(grid.VI(:,1));
      %  a22 = grid.Y(grid.VI(:,3))-grid.Y(grid.VI(:,1));
      %  grid.A = 0.5 * (a11.*a22 - a21.*a12);
        grid.Ainv = grid.A.^(-1);

        % midpoint coordinates of grid-cells
        % equals cog of corners

        % coordinates of vertices: XX(elnum, vnum)
        XX = reshape(grid.X(grid.VI(:)),size(grid.VI));
        YY = reshape(grid.Y(grid.VI(:)),size(grid.VI));
        grid.CX = mean(XX,2);
        grid.CY = mean(YY,2);

        % edge-vectors: [DX(el,edge) DX (el,edge)]
        DX = [XX(:,2:3), XX(:,1)] - XX;
        DY = [YY(:,2:3), YY(:,1)] - YY;

        % matrix with edge-lengths
        grid.EL = sqrt(DX.^2 + DY.^2);

        % matrix with (unit) normal components
        grid.NX =   grid.EL.^(-1) .* DY;
        grid.NY = - grid.EL.^(-1) .* DX;

        % matrix with edge-midpoint-coordinates
        % this computation yields epsilon-differing edge-midpoints on
        % neighbouring elements, perhaps solve later.
        grid.ECX = XX + DX * 0.5;
        grid.ECY = YY + DY * 0.5;

        %%% determine indices of neighbouring elements, default: -1
        NBI = -1 * ones(grid.nelements, grid.nneigh);
        INB = zeros(grid.nelements, grid.nneigh);
        % and now it gets complicated...

        % setup edge list e(1,i), e(2,i) are the point indices of the i-th edge
        %                 eid(i) is the element number, which has the edge
        %                 led(i) : local edge number, which is the edge
        %                          e(:,i) in eid(i)
        ee = zeros(2,grid.nelements*3);
        % li = [1 2 2 3 3 1 4 5 5 6 6 4 ... ]
        li = 1:grid.nelements*3;
        li = [li;li];
        li = reshape(li,6,grid.nelements);
        li = [li(2:end,:); li(1,:)];
        li = reshape(li,2,grid.nelements*3);
        ee = t(li);
        ee = sort(ee); % now smallest vindex above
        elid = repmat(1:grid.nelements,3,1);
        elid = elid(:)';
        led = repmat((1:3)',1,grid.nelements);
        led = led(:)';

        % generate sort key from both indices
        key = ee(1,:) + grid.nvertices * ee(2,:);
        [keysort, ind] = sort(key);
        ee_sort = ee(:,ind);
        elid_sort = elid(ind);
        led_sort = led(ind);

        % find indices of duplicates by simple forward difference
        %  => these are inner edges!

        diff = keysort(1:end-1) - [keysort(2:end)];
        inner_edges = find(diff==0);

        % prepare entries in NBI: for all i in inner_edges
        % NBI(elid_sort(i),led_sort(i)) = elid_sort(i+1);
        % INB(elid_sort(i),led_sort(i)) = led_sort(i+1);
        % und umgekehrt!

        li = sub2ind(size(NBI),...
      	       elid_sort(inner_edges),...
      	       led_sort(inner_edges));
        NBI(li) = elid_sort(inner_edges+1);
        INB(li) = led_sort(inner_edges+1);
        li = sub2ind(size(NBI),...
      	       elid_sort(inner_edges+1),...
      	       led_sort(inner_edges+1));
        NBI(li) = elid_sort(inner_edges);
        INB(li) = led_sort(inner_edges);

      %  % find all element-pairs, which have a common edge
      %  % i.e. neighs(i,1) and neighs(i,2) are two elements, which
      %  % have the same edge i, where this edge is the
      %  % eind(i,1)-th edge (1..3) in the first element
      %  % eind(i,2)-th edge (1..3) in the first element
      %
      %  nedges = max(t(:));
      %  ehist = hist(t(:),1:nedges);
      %  if max(ehist)>3
      %    error('edge occuring in more than 2 triangles!!');
      %  end;
      %  inner_edges = find(ehist==2);
      %  neighs = zeros(length(inner_edges),2);
      %  eind = zeros(length(inner_edges),2);
      %
      %  % somehow map edge-numbers => triangle numbers
      %  % need two of them, as occasionally global and local edge numbers
      %  % can coincide:
      %  % edgemap1(i,locedgenum) = elementno    if edge i is local edge in elementno
      %  % edgemap2(i,locedgenum) = elementno    if edge i is local edge in elementno
      %
      %  edgemap1 = nan * ones(nedges,grid.nneigh);%
      %
      %  [locedgenum,triangnum] = ind2sub(size(t),1:length(t(:)));
      %  li = sub2ind(size(edgemap1),t(:),locedgenum(:));
      %  edgemap1(li) = triangnum; % here due to duplicates, only the last
      %                            % entry is taken
      %
      %  % remove these entries from t and repeat for 2nd edgemap
      %  t2 = t;
      %  [edgenum,locedgenum] = find(~isnan(edgemap1));
      %  li = sub2ind(size(edgemap1),edgenum(:),locedgenum(:));
      %
      %  li2 = sub2ind(size(t2),locedgenum(:),edgemap1(li));
      %  t2(li2) = 0;
      %
      %
      %  edgemap2 = nan * ones(nedges,grid.nneigh);
      %  nonzeros = find(t2~=0);
      %
      %  [locedgenum,triangnum] = ind2sub(size(t2),1:length(nonzeros));
      %  li = sub2ind(size(edgemap2),t2(nonzeros),locedgenum(:));
      %  edgemap2(li) = triangnum;
      %  % now together edgemap1 and edgemap 2 represent the edge information.
      %
      %  % extract inner edges
      %  edgemap = [edgemap1(inner_edges,:), edgemap2(inner_edges,:)];
      %  % check, that exactly two entries per row are non-nan
      %  entries = sum(~isnan(edgemap),2);
      %  if (min(entries)~=2) | max(entries)~=2
      %    disp('error in neighbour detection, please check');
      %    keyboard;
      %  end;
      %
      %  [mi, imi] = min(edgemap,[],2);
      %  [ma, ima] = max(edgemap,[],2);
      %  neighs = [mi, ma];
      %  eind = [imi, ima];
      %
      %  % then fill NBI and INB in a vectorized way, i.e.
      %  % for all edges:
      %  % NBI(neighs(i,1), eind(i,1)) = neighs(i,2);
      %  % NBI(neighs(i,2), eind(i,2)) = neighs(i,1);
      %  % INB(neighs(i,1), eind(i,1)) = eind(i,2);
      %  % INB(neighs(i,2), eind(i,2)) = eind(i,1);
      %
      %  l1 = sub2ind( size(NBI), neighs(:,1) , eind(:,1));
      %  l2 = sub2ind( size(NBI), neighs(:,2) , eind(:,2));
      %  NBI(l1) = neighs(:,2);
      %  NBI(l2) = neighs(:,1);
      %  INB(l1) = eind(:,2);
      %  INB(l2) = eind(:,1);

        % search non-inner boundaries and their midpoints

        li = find(NBI==-1);
        SX = grid.ECX(li(:));
        SY = grid.ECY(li(:));

        % default-boundary = -1
        bnd_ind = -1 * ones(length(li),1);

        if isfield(params,'bnd_rect_index')
          if isfield(params,'boundary_type')
            error(['Do not specify both bnd_rect_index and boundary_type', ...
                   'Äin triagrid construction!']);
          end;
          if (max(params.bnd_rect_index)>0)
            error('boundary indices must be negative!');
          end;
          if size(params.bnd_rect_corner1,1) == 1
            params.bnd_rect_corner1 = params.bnd_rect_corner1';
          end;
          if size(params.bnd_rect_corner2,1) == 1
            params.bnd_rect_corner2 = params.bnd_rect_corner2';
          end;
          for i = 1:length(params.bnd_rect_index)
            indx = (SX > params.bnd_rect_corner1(1,i)) & ...
             (SX < params.bnd_rect_corner2(1,i)) & ...
             (SY > params.bnd_rect_corner1(2,i)) & ...
             (SY < params.bnd_rect_corner2(2,i));
            bnd_ind(indx) = params.bnd_rect_index(i);
          end;
        else
          if isfield(params,'boundary_type')
            bnd_ind = params.boundary_type([SX(:),SY(:)],params);
          end;
        end;

        NBI(li) = bnd_ind'; % set neighbours to boundary

        grid.NBI = NBI;
        grid.INB = INB;

        % check grid consistency:
        nonzero = find(NBI>0); % vector with vector-indices
        [i,j] = ind2sub(size(NBI), nonzero); % vectors with matrix indices
        NBIind = NBI(nonzero); % vector with global neighbour indices
        INBind = INB(nonzero);
        i2 = sub2ind(size(NBI),NBIind, INBind);
        i3 = NBI(i2);
        if ~isequal(i3,i)
      %    plot_element_data(grid,grid.NBI,params);
          disp('neighbour indices are not consistent!!');
          keyboard;
        end;

        % matrix with centroid-distances
        grid.DC = nan * ones(size(grid.CX,1),grid.nneigh);
        nonzero = find(grid.NBI>0);
        [elind, nbind] = ind2sub(size(grid.NBI),nonzero);

        CXX = repmat(grid.CX(:),1,grid.nneigh);
        CYY = repmat(grid.CY(:),1,grid.nneigh);

        CXN = CXX;
        CXN(nonzero) = grid.CX(grid.NBI(nonzero));

        CYN = CYY;
        CYN(nonzero) = grid.CY(grid.NBI(nonzero));

        DC = sqrt((CXX-CXN).^2 + (CYY - CYN).^2);

        grid.DC(nonzero) = DC(nonzero);

        % computation of boundary distances: twice the distance to the border

        nondef = find(isnan(grid.DC));
        nondef = nondef(:)';
        [elind, nind] = ind2sub(size(grid.DC),nondef);
        nind_plus_one = nind + 1;
        i = find(nind_plus_one > grid.nneigh);
        nind_plus_one(i) = 1;
        nondef_plus_one = sub2ind(size(grid.DC),elind,nind_plus_one);
        nondef_plus_one = nondef_plus_one(:)';

        d = zeros(length(nondef),1);

        q = [grid.CX(elind)'; grid.CY(elind)'];
        p1 = [grid.X(grid.VI(nondef)), ...
      	grid.Y(grid.VI(nondef))];
        p2 = [grid.X(grid.VI(nondef_plus_one)), ...
      	grid.Y(grid.VI(nondef_plus_one))];

        grid.DC(nondef) = 2 * dist_point_line(q,p1,p2);

        % make entries of ECX, ECY exactly identical for neighbouring elements!
        % currently by construction a small eps deviation is possible.
        %averaging over all pairs is required
        nonzero = find(grid.NBI>0); % vector with vector-indices
        [i,j] = ind2sub(size(grid.NBI), nonzero); % vectors with matrix indices
        NBIind = NBI(nonzero); % vector with global neighbour indices
        INBind = INB(nonzero); % vector with local edge indices
        i2 = sub2ind(size(NBI),NBIind, INBind);
        % determine maximum difference in edge-midpoints, but exclude
        % symmetry boundaries by relative error < 0.0001
        diffx = abs(grid.ECX(nonzero)-grid.ECX(i2));
        diffy = abs(grid.ECY(nonzero)-grid.ECY(i2));
        fi = find ( (diffx/(max(grid.X)-min(grid.X)) < 0.0001) &  ...
      	      (diffy/(max(grid.Y)-min(grid.Y)) < 0.0001) );

        %disp(max(diffx));
        %disp(max(diffy));
        %  => 0 ! :-)
        % keyboard;

        grid.ECX(nonzero(fi)) = 0.5*(grid.ECX(nonzero(fi))+ grid.ECX(i2(fi)));
        grid.ECY(nonzero(fi)) = 0.5*(grid.ECY(nonzero(fi))+ grid.ECY(i2(fi)));

        % for diffusion discretization: Assumption of points with
        % orthogonal connections to edges. Distances and intersections
        % determined here. For triangles, this can simply be the
        % circumcenters and the corresponding intersections

        q = [grid.X(grid.VI(:,1)), ...
             grid.Y(grid.VI(:,1))];
        p1 = [grid.X(grid.VI(:,2)), ...
              grid.Y(grid.VI(:,2))];
        p2 = [grid.X(grid.VI(:,3)), ...
              grid.Y(grid.VI(:,3))];
        S = circumcenter_triangle(q,p1,p2);
        grid.SX = S(:,1);
        grid.SY = S(:,2);

        % compute DS for inner edges (identical as DC computation!)
        grid.DS = nan * ones(size(grid.CX,1),grid.nneigh);
        nonzero = find(grid.NBI>0);
        [elind, nbind] = ind2sub(size(grid.NBI),nonzero);

        SXX = repmat(grid.SX(:),1,grid.nneigh);
        SYY = repmat(grid.SY(:),1,grid.nneigh);

        SXN = SXX;
        SXN(nonzero) = grid.SX(grid.NBI(nonzero));

        SYN = SYY;
        SYN(nonzero) = grid.SY(grid.NBI(nonzero));

        DS = sqrt((SXX-SXN).^2 + (SYY - SYN).^2);

        grid.DS(nonzero) = DS(nonzero);

        % computation of boundary distances: twice the distance to the border

        nondef = find(isnan(grid.DS));
        nondef = nondef(:)';
        [elind, nind] = ind2sub(size(grid.DS),nondef);
        nind_plus_one = nind + 1;
        i = find(nind_plus_one > grid.nneigh);
        nind_plus_one(i) = 1;
        nondef_plus_one = sub2ind(size(grid.DS),elind,nind_plus_one);
        nondef_plus_one = nondef_plus_one(:)';

        d = zeros(length(nondef),1);

        q = [grid.SX(elind)'; grid.SY(elind)'];
        p1 = [grid.X(grid.VI(nondef)), ...
      	grid.Y(grid.VI(nondef))];
        p2 = [grid.X(grid.VI(nondef_plus_one)), ...
      	grid.Y(grid.VI(nondef_plus_one))];

%        keyboard;

        grid.DS(nondef) = 2 * dist_point_line(q,p1,p2);

	if ~isempty(find(isnan(grid.DS)))
	  error('nans produced in grid generation!');
	end;
	
        % intersections of S_i are the centroids of edges
        grid.ESX = grid.ECX;
        grid.ESY = grid.ECY;

        % compute geometry constants for CFL-condition
        h = max(grid.EL,[],2); % elementwise maximum edgelength
        grid.hmin = min(h);
        % alpha1 * h_j^2 <= |T_j|
        alpha1 = min(grid.A .* h.^(-2));
        % alpha2 * |boundary T_j| <= h_j
        b = sum(grid.EL,2);
        alpha2 = min(h .* b.^(-1));
        % alpha3 h_j <= d_jl, d_jl distance of circumcenters
        alpha3 = min(min(grid.DS,[],2) .* h.^(-1));
        grid.alpha = min([alpha1,alpha2,alpha3]);	

%        if grid.alpha<=0
%          disp(['Caution: Grid-geometry constant alpha less/equal zero', ...
%          ' problematic ' ...
%          ' in automatic CFL-conditions!']);
%        end;

        % DF = [(p2-p1) (p3-p1)] a 2x2 matrix
        DF = zeros(grid.nelements,2,2); % jacobian
        DF(:,1,1) = XX(:,2)-XX(:,1);
        DF(:,2,1) = YY(:,2)-YY(:,1);
        DF(:,1,2) = XX(:,3)-XX(:,1);
        DF(:,2,2) = YY(:,3)-YY(:,1);
        DetDF = DF(:,1,1).*DF(:,2,2)- DF(:,2,1).*DF(:,1,2); % determinant
        DetDFinv = DetDF.^(-1);
        JIT = zeros(grid.nelements,2,2);
        % inversetransposed of A = (a,b; c,d) is A^-1 = 1/det A (d,-c; -b,a);
        JIT(:,1,1)=DF(:,2,2).*DetDFinv;
        JIT(:,1,2)=-DF(:,2,1).*DetDFinv;
        JIT(:,2,1)=-DF(:,1,2).*DetDFinv;
        JIT(:,2,2)=DF(:,1,1).*DetDFinv;

        % test:
        if (size(JIT,1)>=5) &(max(max(abs(...
            transpose(reshape(JIT(5,:,:),2,2))*...
            reshape(DF(5,:,:),2,2)-eye(2)...
            )))>1e-5)
          error('error in JIT computation!');
        end;

        % JIT:
        grid.JIT = JIT; % perhaps later also store area, etc.

        grid.nedges_boundary = length(find(grid.NBI<0));

        grid.nedges_interior = 0.5*(3*grid.nelements - grid.nedges_boundary);
        if ~isequal(ceil(grid.nedges_interior), ...
                    grid.nedges_interior);
          error('number of inner edges odd!!!');
        end;

      end
    end

    demo(dummy);

    display(grid);

    gridp = gridpart(grid,eind);

    lcoord = llocal2local(grid, faceinds, llcoord);

    glob = local2global(grid, einds, loc, params);

    p = plot(gird, params);

    grid = set_nbi(grid, nbind, values);

    function gcopy = copy(grid)
      % @copybrief gridbase::copy
      %
      % Return values:
      %   gcopy: object of type triagrid. This is a deep copy of the current instance.
      %

      gcopy = triagrid(grid);
    end

  end

  methods(Static)

      [C, G] = aff_trafo_glob2loc(x0, y0);

      [C, G] = aff_trafo_loc2glob(x0, y0);

      [C, G] = aff_trafo_orig2ref(x0, y0, varargin);

      loc = global2local(grid, elementid, glob);

      micro2macro = micro2macro_map(microgrid, macrogrid);
  end
end

