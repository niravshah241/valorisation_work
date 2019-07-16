function [res] = pressure_average_velocity_basis_jump_int_minus_int_plus...
    (grid,params,paramsP,tria_index,local_vertex_index,qdeg)

f = @(llcoord) pressure_average_velocity_basis_jump_local_minus_local_plus...
    (llcoord, grid,params,paramsP, tria_index,local_vertex_index);

face_index = local_vertex_index;

res = intervalquadrature(qdeg,f) * grid.EL(tria_index,face_index);