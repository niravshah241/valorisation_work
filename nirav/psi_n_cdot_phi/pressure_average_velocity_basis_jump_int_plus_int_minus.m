function [res] = pressure_average_velocity_basis_jump_int_plus_int_minus...
    (grid,params,paramsP,tria_index,local_vertex_index,qdeg)

f = @(llcoord) pressure_average_velocity_basis_jump_local_plus_local_minus...
    (llcoord, grid,params,paramsP, tria_index,local_vertex_index);

face_index = local_vertex_index;

res = intervalquadrature(qdeg,f) * grid.EL(tria_index,face_index);