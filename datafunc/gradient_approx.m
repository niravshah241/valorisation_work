function grad = gradient_approx(model,model_data,U, NU_ind, edge)
%function grad = gradient_approx(model,model_data,U, [NU_ind], edge)
%
% function that computes an approximation of the gradient over the edge
% given by 'edge'.
%
% Martin Drohmann 13.05.2008

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


  grid = model_data.grid;
  nels = grid.nelements;

  if(isempty(NU_ind))
    NU_ind = 1:nels;
  end
  % get a 6-environment of each cell (see get_enbi for details)
  ENBI    = get_enbi(grid, edge, model.tstep);
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
  if(any(any(ENBI(NU_ind, :) == -10)))
    disp('Attention in gradient_approx! This line should never be executed');
  end
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
%  assert(max(max(ENBI <= 0)) == 0);

  UENBI = U(ENBI);
  % + ytrans(grid.CX(ENBI(non_zero_ind)), grid.CY(ENBI(non_zero_ind))); 

  % construct the 4-environment UE            
  UE(NU_ind,[1,2]) = UENBI(NU_ind,[1,4]);
  UE(NU_ind,3)     = sum(UENBI(NU_ind,[1,2,4,5]), 2) / 4;
  UE(NU_ind,4)     = sum(UENBI(NU_ind,[1,3,4,6]), 2) / 4;
  % update the 4-environment concerning dirichlet cells
  Xdir = grid.ESX(case1_dir_ind, edge);
  Ydir = grid.ESY(case1_dir_ind, edge);
  UE(case1_dir_ind, 2) = model.dirichlet_values_ptr([Xdir,Ydir], model);
  Xdir = grid.X(grid.VI(case2_dir_ind, case23_mapper(1) ));
  Ydir = grid.Y(grid.VI(case2_dir_ind, case23_mapper(1) ));
  UE(case2_dir_ind, 3) = model.dirichlet_values_ptr([Xdir,Ydir], model);
  Xdir = grid.X(grid.VI(case3_dir_ind, case23_mapper(2) ));
  Ydir = grid.Y(grid.VI(case3_dir_ind, case23_mapper(2) ));
  UE(case3_dir_ind, 4) = model.dirichlet_values_ptr([Xdir,Ydir], model);

  horiz_length = grid.EL(1,1);
  vert_length  = grid.EL(1,2);

  if edge == 1
    grad = [(UE(NU_ind,2)-UE(NU_ind,1))./horiz_length, (UE(NU_ind,4)-UE(NU_ind,3))./vert_length];
  elseif edge == 2
    grad = [(UE(NU_ind,3)-UE(NU_ind,4))./horiz_length,  (UE(NU_ind,2)-UE(NU_ind,1))./vert_length];
  elseif edge == 3
    grad = [(UE(NU_ind,1)-UE(NU_ind,2))./horiz_length, (UE(NU_ind,4)-UE(NU_ind,3))./vert_length];
  elseif edge == 4
    grad = [(UE(NU_ind,3)-UE(NU_ind,4))./horiz_length,  (UE(NU_ind,1)-UE(NU_ind,2))./vert_length];
  end

%| \docupdate 
