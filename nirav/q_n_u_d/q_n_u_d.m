function [res]=q_n_u_d(lcoord,params,paramsP,grid, tria_index, local_vertex_index)
% L^2 scalar product (q,n \cdot u_D)
% grid needs fields for normal(grid.NX and grid.NY), local2global, llocal2local and params,paramsP needs
% fields for ldg_evaluate_basis, local2global, llocal2local and tria_index, 
% local_vertex_index are dirichlet element number and dirichlet edge number

%face_index = local_vertex_index;
N=[grid.NX(tria_index,local_vertex_index) grid.NY(tria_index,local_vertex_index)];
%lcoord = llocal2local(grid,local_vertex_index,llcoord);
q=ldg_evaluate_basis(lcoord,paramsP);
u_d=dirichlet_values(local2global(grid,tria_index,lcoord,paramsP),params);
u_d = N*u_d;

res = zeros(paramsP.ndofs_per_element,1);

res=q*u_d;

end