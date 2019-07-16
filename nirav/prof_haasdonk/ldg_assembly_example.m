function A = ldg_assembly_example

% function demonstrating the assembly of a mass matrix for a

% vectorial or scalar ldg function



% B. Haasdonk 11.8.2017



qdeg = 4;

% quadrature degree for integration:

% initialize grid

grid = triagrid();



disp('');

disp('initialization of zero ldg function');



% initialize basic data



params = [];

params.pdeg = 1;

params.dimrange = 2; % for vectorial

%params.dimrange = 1; % for scalar

params.qdeg = qdeg;

% vectorial function



df_info = ldginfo(params, grid);

%params.nelements = grid.nelements;

%params.ndofs = ldg_ndofs(params);

%params.ndofs_per_element = ldg_ndofs_per_element(params);

dofs = ldg_zero(df_info); % simple zero vector of dofs



df = ldgdiscfunc(dofs,df_info); % setting of data: persistent!



% example of mass matrix of ldg space, simple block-diagonal matrix

M = ldg_mass_matrix(df,grid,params);

spy(M);



disp('press enter to continue')



% example of weighted mass matrix of type

%

%  (A)_ij =  int_Omega   a(x) (phi_i(x) * phi_j(x)) dx



% accept glob to be a npoints x 2 matrix, then return column vector [a(glob(1,:)),...a(glob(end,:))]:

a = @(glob,params) ones(size(glob,1),1); % simple constant 1 function



% For the integration need to access to a via elements and local coordinates

a_local = @(einds,loc,grid,params) a(...
    local2global(grid,einds,loc,params),params);



% Denote elements of Mesh by e_k, k=1,...,nelements coordinate

% variable x

% Denote reference triangle as T with coordinate variable y

% and reference map F_k : T -> e_k, F(y) = x

% and Jacobian J_k := D F_k

% and determinant det(J_k)

% Then

%  (A)_ij =  sum_k int_{e_k} a(x)  phi_i(x) * phi_j(x) dx

%         =  sum_k int_{T}   a(F_k(y)) phi_i(F_k(y)) * phi_j(F_k(y)) |det J_k| dy

%         =  sum_k int_{T}   a(F_k(y)) phihat_ihat(y) * phihat_jhat(y) |det J_k| dy

%

%   with ihat and jhat the local dof index of the global functions

%   with number i and j.

% The assembly now is not done by double loop over ij, but single

% loop over the elements, computing the local element contributions

%   A_local   =  int_{T} a(F_k(y)) phihat_ihat(y) * phihat_jhat(y) |det J_k| dy

% and then adding this matrix to the global matrix A:



A = sparse(df.ndofs, df.ndofs);



gids = ldg_global_dof_index(df,grid);

% accept single loop over mesh (vectorize sometime later)

for k = 1:grid.nelements;
    
    abs_det_J_k = grid.A(k)*2; % determininat is parallelogram hence twice triangle
    
    
    
    % assemble local weighted mass matrix
    
    f = @(lcoord) a_local(k,lcoord,grid,params) * gram_matrix(ldg_evaluate_basis(lcoord,df_info)');
    
    A_local = triaquadrature(qdeg,f) * abs_det_J_k;
    
    
    
    % distribute local matrix into global matrix by dofmap
    
    lids = 1:df.ndofs_per_element; % local dof indices
    
    elid = k; % element index
    
    ids = gids(k,:);
    
    A(ids,ids) = A(ids,ids) + A_local;
    
end;



% check that for a(x)==1 the above resembles the mass matrix M



err = max(max(abs(A - M)));

if err<1e-10
    
    disp(['Test OK: assembly is corresponding with simple block-diagonal mass' ...
    ' matrix']);

else
    
    disp('Error in assembly!!!');
    
end;





