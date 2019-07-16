function grad = gradient_approx2(U, grid, params, edge)

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

  nels = grid.nelements;


  % get a 6-environment of each cell (c.f. get_ENBI for details)
  ENBI    = get_ENBI(grid, edge);
  % initialise the U values of the 6 environment
%  UENBI   = zeros(size(ENBI));
  % initialise UE, a 4 environment that stores interpolated values at the
  % edges tips and exact values at the two adjacent cell centers. From UE
  % the gradient is calculated by lineare interpolation (lsq)
  UE          = zeros(nels, 4);

  %%%% Handle dirichlet cells ("ghost" cells with index -1) %%%%
  % if a cell adjacent to the current edge is a _dirichlet_ "ghost" cell,
  % entries in the matrix UE have to be replaced by exact evaluations at
  % the dirichlet boundary. The following cases have to be separated:
  % case 1: one of the adjacent cells is a ghost
  %      => replace the cell's entry in UE with the dirichlet value at
  %         the edge's center
  % case 2/3: one cell adjacent to the tip is a ghost
  %      => replace the tip's entry in UE with the dirichlet value at the
  %         edge's tip.

  % find the _row_ index (edge inedex) of the three cases
  case1_dir_ind = find(ENBI(:,4) == -1);
  case2_dir_ind = find(ENBI(:,5) == -1);
  case3_dir_ind = find(ENBI(:,6) == -1);

  % the case23_mapper maps to the correct edge tip (correct column index)
  % in UE.
  if edge > 2
    case23_mapper = [ mod(edge, 4) + 1, edge ];
  else
    case23_mapper = [ edge, mod(edge, 4) + 1 ];
  end

  % neuman boundaries!!! For all other ghost cells the entries of
  % adjacent non-ghost cells are copied. Ghost cells in the corners of
  % the rectgrid are assigned in the end by copying one of the already
  % assigned adjacent ghost-cell. It is unspecified which ghost cell is
  % used.
  nb_ind          = find(ENBI(:,2) <= 0);
  ENBI(nb_ind, 2) = ENBI(nb_ind, 1);
  nb_ind          = find(ENBI(:,3) <= 0);
  ENBI(nb_ind, 3) = ENBI(nb_ind, 1);
  nb_ind          = find(ENBI(:,4) <= 0);
  ENBI(nb_ind, 4) = ENBI(nb_ind, 1);
  nb_ind          = find(ENBI(:,5) == -2);
  ENBI(nb_ind, 5) = ENBI(nb_ind, 2);
  nb_ind          = find(ENBI(:,5) <= 0);
  ENBI(nb_ind, 5) = ENBI(nb_ind, 4);
  nb_ind          = find(ENBI(:,6) == -2);
  ENBI(nb_ind, 6) = ENBI(nb_ind, 3);
  nb_ind          = find(ENBI(:,6) <= 0);
  ENBI(nb_ind, 6) = ENBI(nb_ind, 4);

  % now every ghost cell should have a positive index
  %assert(max(max(ENBI <= 0)) == 0);

  UENBI = U(ENBI);
  % + ytrans(grid.CX(ENBI(non_zero_ind)), grid.CY(ENBI(non_zero_ind))); 

  % construct the 4-environment UE            
  UE(:,[1,2]) = UENBI(:,[1,4]);
  UE(:,3)     = sum(UENBI(:,[1,2,4,5]), 2) / 4;
  UE(:,4)     = sum(UENBI(:,[1,3,4,6]), 2) / 4;

  % construct the 4-environment UE            
  UE(:,[1,2]) = UENBI(:,[1,4]);
  UE(:,3)     = sum(UENBI(:,[1,2,4,5]), 2) / 4;
  UE(:,4)     = sum(UENBI(:,[1,3,4,6]), 2) / 4;
  % update the 4-environment concerning dirichlet cells
  Xdir = grid.ESX(case1_dir_ind, edge);
  Ydir = grid.ESY(case1_dir_ind, edge);
  %        Yoff = ytrans(Xdir, Ydir);
  UE(case1_dir_ind, 2) = dirichlet_values(Xdir, Ydir, params);% + Yoff;

  Xdir = grid.X(grid.VI(case2_dir_ind, case23_mapper(1) ));
  Ydir = grid.Y(grid.VI(case2_dir_ind, case23_mapper(1) ));
  %        Yoff = ytrans(Xdir, Ydir);
  UE(case2_dir_ind, 3) = dirichlet_values(Xdir, Ydir, params);% + Yoff;
  Xdir = grid.X(grid.VI(case3_dir_ind, case23_mapper(2) ));
  Ydir = grid.Y(grid.VI(case3_dir_ind, case23_mapper(2) ));
  %        Yoff = ytrans(Xdir, Ydir);
  UE(case3_dir_ind, 4) = dirichlet_values(Xdir, Ydir, params);% + Yoff;

  % construct the LHS matrix of the LES, it is - up to the sign -
  % identical for edges 1 and 3, respectively 2 and 4.
  % Therefore we need this factor 'factor' and an dummy edge 'tmp':
  factor = -1;
  tmp    = edge;
  if(edge > 2)
    tmp    = edge - 2;
    factor = 1;
  end

  nb_ind = grid.NBI(1, tmp);
  Acell  = [ grid.CX([1; nb_ind])           , grid.CY([1; nb_ind]); ...
             grid.X(grid.VI(1,(0:1) + tmp))', grid.Y(grid.VI(1,(0:1) + tmp))'];
  Acell  = ( Acell - repmat([grid.ECX(1,tmp), grid.ECY(1, tmp)], 4, 1) ) ...
             * factor;

  AcellT   = Acell';
  lhs_cell = AcellT * Acell;
  % The 3 diagonals of the LHS matrix A^T*A.
  LHS_d    = repmat(horzcat([ lhs_cell(2,1); 0 ], ...
                    diag(lhs_cell), ...
                    [ 0; lhs_cell(1,2) ]),...
                    nels, 1);
  LHS      = spdiags(LHS_d, -1:1, 2*nels, 2*nels);

  % The sparse matrix A^T (consists of 2x4 matrix blocks)
  AT       = sparse( repmat([1;2], 4*nels, 1) ...
                     + reshape(repmat(0:2:2*nels-1, 8, 1), 8*nels, 1), ...
                     reshape(repmat(1:4*nels, 2, 1), 8*nels, 1), ...
                     repmat(reshape(AcellT, 8, 1), nels, 1) );
%  full(AT)
  % RHS = A^T b                   
  RHS      = AT * reshape(UE', 4*nels, 1);

  grad  = reshape(LHS \ RHS, 2, nels)';
  %grad
  %keyboard


% TO BE ADJUSTED TO NEW SYNTAX
%| \docupdate 
