function U = fem_laplace(grid,params)
%function U = fem_laplace(grid,params)
% 
% function solving the elliptic laplace problem with general dirichlet and
% neuman boundary data using a nodal basis. currently cartesian
% grid (rectgrid) is assumed and bilinear lagrange-elements. More general
% ansatz-functions and quadratures should be implemented in the future.
%
% Current problem: discrete gradient of the solution is not divergence
% free so must be divergence-cleaned if used as velocity field 
% in other simulations

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

  
% Bernard Haasdonk 13.4.2006:

%if isempty(grid)
%  grid=grid_geometry(params);
%end;

if ~isa(grid,'rectgrid')
  error('currently only rectgrid supported for fem!!');
end;  

% indices of points on dirichlet boundary  
  dir_vertex = zeros(size(grid.X));

  if params.verbose > 5
    disp('entering fe_laplace');
  end;
  
  % search all dirichlet edges and mark both corresponding vertices
  if params.verbose > 5
    disp('searching dirichlet edges');
  end;
  dir_edge = grid.NBI == -1;
  nfaces = size(grid.NBI,2);
  for i = 1:nfaces;
    di = find(dir_edge(:,i));
    dir_vertex(grid.VI(di,i)) = 1;
    dir_vertex(grid.VI(di,mod(i,nfaces)+1)) = 1;
  end;
  dir_vind = find(dir_vertex);
  
  ndof = length(grid.X);
  dirichlet_fct = zeros(ndof,1);
  
  % determine dirichlet boundary values
  dirichlet_fct(dir_vind) = dirichlet_values(grid.X(dir_vind),...
					     grid.Y(dir_vind),params);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % setup stiffness matrix as a sparse matrix
  % S(i,j) = 1 if i==j and x_i is a dirichlet-vertex
  % S(i,j) = 0 if i!=j and x_i is a dirichlet-vertex
  % S(i,j) = int sigma * grad phi_i * grad_phi_j if x_i,j are non-dirichlet  
  %
  % so first compute complete stiffness-matrix
  % (then kronecker-kill later on before system solving)
  if params.verbose > 5
    disp('setting up stiffness matrix');
  end;
  %S = sparse([0]zeros(ndof));
  S = spalloc(ndof,ndof,13*ndof);
  
  % in the general case: compute all local matrices simultaneously
  % L(i,j,k) contribution of element i, local basis functions j and k 
  %
  % in our case, L is not element-dependent, so only compute L(j,k)
  %
  % for cartesian grid and constant sigma the local matrix is
  % identical for all elements
  % domain: [0, dx] x [0, dy], 
  % counting right lower corner 1 then counterclockwise
  % detAinv = (dx*dy)^-1
  % grad phi_1 = detAinv * [dy-y, -x   ]' 
  % grad phi_2 = detAinv * [   y,  x   ]' 
  % grad phi_3 = detAinv * [  -y,  dx-x]' 
  % grad phi_4 = detAinv * [y-dy,  x-dx]' 
  % =>  L(j,k) = detAinv * sigma
  %   [1/3dx^2+1/3dy^2)  -1/3dx^2+1/6dy^2  -1/6dx^2-1/6dy^2   1/6dx^2-1/3dy^2 
  %    XXXX               1/3dx^2+1/3dy^2   1/6dx^2-1/3dy^2  -1/6dx^2-1/6dy^2 
  %    XXXX                 XXXX            1/3dx^2+1/3dy^2  -1/3dx^2+1/6dy^2 
  %    xxxx                 XXXX             XXXX             1/3(dx^2+dy^2) ]
  % 
  nbasisfcts = size(grid.VI,2);
  sigma = 1;
  dx = (params.xrange(2)-params.xrange(1))/params.xnumintervals;
  dy = (params.yrange(2)-params.yrange(1))/params.ynumintervals;
  
  Lx = [1/3 -1/3 -1/6  1/6; 
       -1/3  1/3  1/6 -1/6;
       -1/6  1/6  1/3 -1/3;
        1/6 -1/6 -1/3  1/3];

  Ly = [1/3  1/6  -1/6 -1/3;
        1/6  1/3  -1/3 -1/6;
       -1/6 -1/3   1/3  1/6;
       -1/3 -1/6   1/6  1/3];

  L = (dx*dy)^(-1) * sigma * (Lx * dx^2 + Ly * dy^2);
      
  % accumulate local matrices to single operator matrix
  % avoid large element loop, instead perform loop over pairs of
  % local basis functions
  
  for j= 1:nbasisfcts
    for k= 1:nbasisfcts
      if params.verbose > 5
	disp(['(j,k)=(',num2str(j),',',num2str(k),')']);
      end;
      global_dofs_j = grid.VI(:,j); 
      global_dofs_k = grid.VI(:,k);       
      global_dofs_ind = sub2ind(size(S),global_dofs_j,global_dofs_k);
      if max(global_dofs_ind)>1.6e9
	% the following code results in t = 707.5058 s for 300x150 grid
	%   cut S into column-block-pieces for ommiting large-index accesses
	%	v = lget(S,global_dofs_ind);
	%	v = v + L(j,k);
	%	S = lset(S,global_dofs_ind,v);
	% the following code results in t = 476.4325 s for 300x150 grid
	%   for z = 1:length(global_dofs_ind)
	%     S(global_dofs_j(z),global_dofs_k(z)) =  ...
	%      S(global_dofs_j(z),global_dofs_k(z)) + L(j,k);
	%   end;
	% the following code results in t = 62s s for 300x150 grid :-)
	% most time spent in solution of LGS
	V = sparse(global_dofs_j, global_dofs_k, ...
		   L(j,k),size(S,1),size(S,2));
	S = S + V;
      else % direct subvector operation is possible
	S(global_dofs_ind) = S(global_dofs_ind) + L(j,k); 
	%      S(global_dofs_ind) = S(global_dofs_ind) + L(:,j,k); 	
      end;
    end;
  end;

% alternative: Loop over elements, (which is to be omitted)
%  nelements = length(grid.A);
%  for i = 1:nelements
%  global_dofs_i = ...
%    S(global_dofs_i,global_dofs_i) =  ...
%	S(global_dofs_i,global_dofs_i) + L(i,:,:);   
%  end;
    
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % setup right hand side
  % b(j) = T1 + T2 + T3  if x_j is a non-dirichlet vertex
  % b(j) = 0 if x_j is a dirichlet vertex
  %
  % so first compute all nontrivial quantities, then set dir_vind to 0
  %
  %%%%%%%%% source term:
  %    T1 = int source * phi_j                      (2D integral)
  %              -> zero in our case as we have no source
  if params.verbose > 5
    disp('setting up right hand side');
  end;

  T1 = zeros(ndof,1);

  %%%%%%%%% dirichlet term:
  %    T2 = - int sigma * grad b_dir * grad phi_j   (negative of a 2D integral)
  % after projection of b_dir to the discrete function space 
  % this is exactly the stiffness-matrix multiplied with the b_dir
  % node values
  T2 = - S * dirichlet_fct;

  %%%%%%%%% neuman term:
  %    T3 = int_gamma_neu  b_neu * phi_j            (1D integral)
  %
  % goal: Assume b_neu also piecewise linear -> exact integration possible 
  % for all edges, the entries of the two corresponding nodes are modified
  T3 = zeros(ndof,1);  

  % search all neuman edges and generate list of start and end points
  neu_edge = find(grid.NBI == -2);
  [i,j] = ind2sub(size(grid.NBI),neu_edge);
  % list of global indices of start points
  neu_vind1 = grid.VI(neu_edge);
  % list of global indices of end points
  j2 = mod(j,nfaces)+1; 
  neu_edge_plus1 = sub2ind(size(grid.NBI),i,j2);
  neu_vind2 = grid.VI(neu_edge_plus1);
  
  neu_vals1= neuman_values(grid.X(neu_vind1),...
			   grid.Y(neu_vind1),[],[],[],params);
  neu_vals2= neuman_values(grid.X(neu_vind2),...
			   grid.Y(neu_vind2),[],[],[],params);
  % generate lists of neuman integral values
  % neu_int1(neu_edgenr) is the integral contribution to start node 
  neu_lengths = grid.EL(neu_edge);

  % determine elementary integrals 
  % (phi_1 is basis function of start point, phi_2 of end point) 
  % int_edge Phi_1 Phi_1 = 1/3 |S|
  % int_edge Phi_1 Phi_2 = 1/6 |S|
  % int_edge Phi_2 Phi_2 = 1/3 |S|
  
  % neu_int1 = int_edge Phi_1 (neu_val1 * Phi1 + neu_val2 * Phi2)
  neu_int1 = neu_lengths.* (1/3 * neu_vals1 + 1/6 * neu_vals2); 
  neu_int2 = neu_lengths.* (1/6 * neu_vals1 + 1/3 * neu_vals2); 
  T3(neu_vind1) = T3(neu_vind1) + neu_int1; 
  T3(neu_vind2) = T3(neu_vind2) + neu_int2; 
    
  % assemble
  b = T1 + T2 + T3;  
  b(dir_vind) = 0;
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %delete dirichlet-nodes in matrix
  S(dir_vind,:) = 0;
  S(:,dir_vind) = 0;
  S(dir_vind,dir_vind) = speye(length(dir_vind));
  
  % solve homogeneous problem  
  % U = S^-1 * b;
  if params.verbose > 5
    disp('solving LGS system');
  end;
    
  %[U, flag] = pcg(S,b);
%  keyboard;
% nonsymmetric solvers:
%  [U, flag] = bicgstab(S,b,[],1000); => 65sec
%  [U, flag] = cgs(S,b,[],1000); % => 47 sec
% symmetric solver, non pd:
  [U, flag] = symmlq(S,b,[],1000); % => 31 sec
%  [U, flag] = minres(S,b,[],1000); % => 31 sec
% symmetric solver, pd:
%  [U, flag] = pcg(S,b,[],1000); % => 46 sec
  if flag>0
    disp(['error in system solution, solver return flag = ', ...
	   num2str(flag)]);
    keyboard;
  end;
  
  % add dirichlet boundary values
  U = U + dirichlet_fct;  
  

  
  

% TO BE ADJUSTED TO NEW SYNTAX
%| \docupdate 
