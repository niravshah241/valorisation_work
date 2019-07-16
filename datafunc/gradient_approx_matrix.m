function [grad,bdir] = gradient_approx_matrix(model,model_data, NU_ind, edge)
%function grad = gradient_approx_matrix(model,model_data, [NU_ind], edge)
%
% function that returns a matrix with which an approximation of the gradient
% over the edge given by 'edge' can be computed.
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


decomp_mode = model.decomp_mode;

if decomp_mode == 2
  grad = 1;
  bdir = model.dirichlet_values_ptr([],model);
else
  grid = model_data.grid;
  nels = grid.nelements;

  if(isempty(NU_ind))
    NU_ind = 1:nels;
  end
  nu_ind_length = length(NU_ind);



  % get a 6-environment of each cell (see get_enbi for details)
  ENBI = get_enbi(grid, edge, model.tstep);

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

  % find the _row_ index (edge index) of the three cases
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

  % update the 4-environment concerning dirichlet cells
  if decomp_mode == 0
    UE          = zeros(nels, 3);
    Xdir = grid.ESX(case1_dir_ind, edge);
    Ydir = grid.ESY(case1_dir_ind, edge);
    UE(case1_dir_ind, 1) = model.dirichlet_values_ptr([Xdir,Ydir],model);
    Xdir = grid.X(grid.VI(case2_dir_ind, case23_mapper(1) ));
    Ydir = grid.Y(grid.VI(case2_dir_ind, case23_mapper(1) ));
    UE(case2_dir_ind, 2) = model.dirichlet_values_ptr([Xdir,Ydir],model);
    Xdir = grid.X(grid.VI(case3_dir_ind, case23_mapper(2) ));
    Ydir = grid.Y(grid.VI(case3_dir_ind, case23_mapper(2) ));
    UE(case3_dir_ind, 3) = model.dirichlet_values_ptr([Xdir,Ydir],model);
  else % decomp_mode == 1
    Xdir = grid.ESX(case1_dir_ind, edge);
    Ydir = grid.ESY(case1_dir_ind, edge);
    tmp = model.dirichlet_values_ptr([Xdir,Ydir],model);
    UE = cell(length(tmp), 3);
    UE(:,1) = cellfun(@(X)(assemble_UE(X,nels,case1_dir_ind)), tmp, ...
                      'UniformOutput', false);
    Xdir = grid.X(grid.VI(case2_dir_ind, case23_mapper(1) ));
    Ydir = grid.Y(grid.VI(case2_dir_ind, case23_mapper(1) ));
    tmp = model.dirichlet_values_ptr([Xdir,Ydir],model);
    UE(:,2) = cellfun(@(X)(assemble_UE(X,nels,case2_dir_ind)), tmp, ...
                      'UniformOutput', false);
    Xdir = grid.X(grid.VI(case3_dir_ind, case23_mapper(2) ));
    Ydir = grid.Y(grid.VI(case3_dir_ind, case23_mapper(2) ));
    tmp = model.dirichlet_values_ptr([Xdir,Ydir],model);
    UE(:,3) = cellfun(@(X)(assemble_UE(X,nels,case3_dir_ind)), tmp, ...
                      'UniformOutput', false);
  end

  % neuman boundaries!!! For all other ghost cells the entries of
  % adjacent non-ghost cells are copied. Ghost cells in the corners of
  % the rectgrid are assigned in the end by copying one of the already
  % assigned adjacent ghost-cell. It is unspecified which ghost cell is
  % used.
  if(any(any(ENBI(NU_ind, :) == -10)))
    disp('Attention in gradient_approx_matrix! This line should never be executed');
  end
  dir_bnd_ind1    = ENBI(NU_ind,4) == -1;
  dir_bnd_ind2    = find(ENBI(NU_ind,[2,5]) == -1);
  dir_bnd_ind3    = find(ENBI(NU_ind,[3,6]) == -1);
  case2_dir_ind = find(ENBI(NU_ind,5) == -1);
  case3_dir_ind = find(ENBI(NU_ind,6) == -1);
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

  I = cell(2,1);
  J = cell(2,1);
  S = cell(2,1);

  I{1} = ones(nu_ind_length, 2);
  J{1} = ones(nu_ind_length, 2);
  S{1} = ones(nu_ind_length, 2);
  I{2} = ones(nu_ind_length, 4);
  J{2} = ones(nu_ind_length, 4);
  S{2} = ones(nu_ind_length, 4);
  I{3} = ones(length(case2_dir_ind), 2);
  J{3} = ones(length(case2_dir_ind), 2);
  S{3} = ones(length(case2_dir_ind), 2);
  I{4} = ones(length(case3_dir_ind), 2);
  J{4} = ones(length(case3_dir_ind), 2);
  S{4} = ones(length(case3_dir_ind), 2);
  
  horiz_length = grid.EL(1,1);
  vert_length  = grid.EL(1,2);

  if edge == 1
    s1factor      = -1/horiz_length;
    s2factor      = -1/vert_length;
    order_addend1 = -1;
    order_addend2 = 0;
  elseif edge == 2
    s1factor      = -1/vert_length;
    s2factor      = 1/horiz_length;
    order_addend1 = 0;
    order_addend2 = -1;
  elseif edge == 3
    s1factor      = 1/horiz_length;
    s2factor      = -1/vert_length;
    order_addend1 = -1;
    order_addend2 = 0;
  elseif edge == 4
    s1factor      = 1/vert_length;
    s2factor      = 1/horiz_length;
    order_addend1 = 0;
    order_addend2 = -1;
  end

  I{1} = repmat(2*(1:nu_ind_length)+order_addend1,1,2);
  J{1} = ENBI(NU_ind,[1,4]);
%  J{1}(dir_bnd_ind1,2) = dir_bnd_ind1;
  S{1} = s1factor * repmat([1,-1],nu_ind_length,1);
  S{1}(dir_bnd_ind1,2) = 0;

%  S{1} = 1;
  I{2} = repmat(2*(1:nu_ind_length)+order_addend2,1,4);
  J{2} = ENBI(NU_ind,[2,5,3,6]);
  [row_numbers1,dummy] = ind2sub([nu_ind_length,2],dir_bnd_ind2);
  [row_numbers2,dummy] = ind2sub([nu_ind_length,2],dir_bnd_ind3);
  S{2} = s2factor * repmat([1/4,1/4,-1/4,-1/4],nu_ind_length,1);
  S{2}(row_numbers1,[1,2]) = 0;
  S{2}(row_numbers2,[3,4]) = 0;

  I{3} = repmat(2*case2_dir_ind+order_addend2,1,2);
  J{3} = ENBI(case2_dir_ind,[1,4]);
  S{3} = s2factor * repmat([-1/4,-1/4],length(case2_dir_ind),1);

  I{4} = repmat(2*case3_dir_ind+order_addend2,1,2);
  J{4} = ENBI(case3_dir_ind,[1,4]);
  S{4} = s2factor * repmat([1/4,1/4],length(case3_dir_ind),1);


  II   = [ reshape(I{1},2*nu_ind_length,1); reshape(I{2},4*nu_ind_length,1); ...
           reshape(I{3},2*length(case2_dir_ind),1); reshape(I{4},2*length(case3_dir_ind),1) ];
  JJ   = [ reshape(J{1},2*nu_ind_length,1); reshape(J{2},4*nu_ind_length,1); ...
           reshape(J{3},2*length(case2_dir_ind),1); reshape(J{4},2*length(case3_dir_ind),1) ];
  SS   = [ reshape(S{1},2*nu_ind_length,1); reshape(S{2},4*nu_ind_length,1); ...
           reshape(S{3},2*length(case2_dir_ind),1); reshape(S{4},2*length(case3_dir_ind),1) ];

  if decomp_mode == 0
    if edge == 1
      bdir = [UE(NU_ind,1)./horiz_length, (UE(NU_ind,3)-UE(NU_ind,2))./(vert_length)];
    elseif edge == 2
      bdir = [(UE(NU_ind,2)-UE(NU_ind,3))./(horiz_length),  UE(NU_ind,1)./vert_length];
    elseif edge == 3
      bdir = [-UE(NU_ind,1)./horiz_length, (UE(NU_ind,3)-UE(NU_ind,2))./(vert_length)];
    elseif edge == 4
      bdir = [(UE(NU_ind,2)-UE(NU_ind,3))./(horiz_length),  -UE(NU_ind,1)./vert_length];
    end
    bdir = reshape(bdir', 2*nu_ind_length,1);
    grad = sparse(II,JJ,SS, 2*nu_ind_length,nels);
  else %decomp_mode == 1
    if edge == 1
      tmp1 = cellfun(@(X)(X(NU_ind)./horiz_length), UE(:,1), ...
                     'UniformOutput', false);
      tmp2 = cellfun(@(X,Y)( (X(NU_ind)-Y(NU_ind))./vert_length ), UE(:,3), UE(:,2), ...
                     'UniformOutput', false);
    elseif edge == 2
      tmp1 = cellfun(@(X,Y)( (X(NU_ind)-Y(NU_ind))./horiz_length ), UE(:,2), UE(:,3), ...
                     'UniformOutput', false);
      tmp2 = cellfun(@(X)(X(NU_ind)./vert_length), UE(:,1), ...
                     'UniformOutput', false);
    elseif edge == 3
      tmp1 = cellfun(@(X)(-X(NU_ind)./horiz_length), UE(:,1), ...
                     'UniformOutput', false);
      tmp2 = cellfun(@(X,Y)( (X(NU_ind)-Y(NU_ind))./vert_length ), UE(:,3), UE(:,2),...
                     'UniformOutput', false);
    elseif edge == 4
      tmp1 = cellfun(@(X,Y)( (X(NU_ind)-Y(NU_ind))./horiz_length ), UE(:,2), UE(:,3), ...
                     'UniformOutput', false);
      tmp2 = cellfun(@(X)(-X(NU_ind)./vert_length), UE(:,1), ...
                     'UniformOutput', false);
    end
    bdir = cellfun(@(X,Y)([X,Y]), tmp1, tmp2, 'UniformOutput', false);
    bdir = cellfun(@(X)(reshape(X', 2*nu_ind_length,1)), bdir, 'UniformOutput', false);
    grad = {sparse(II,JJ,SS,2*nu_ind_length,nels)};
  end
end

end

function [ret] = assemble_UE(X,nels,dir_ind)
  ret          = zeros(1,nels);
  ret(dir_ind) = X;
end

%| \docupdate 
