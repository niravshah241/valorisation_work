function [res] = ...
    pressure_average_velocity_basis_jump_int_minus_int_plus_11...
    (grid,params,paramsP,tria_index,local_vertex_index,qdeg)

if nargin == 5
    qdeg = params.qdeg;
end

f = @(llcoord) pressure_average_velocity_basis_jump_local_minus_local_plus_11...
    (llcoord, grid,params,paramsP, tria_index,local_vertex_index);

face_index = local_vertex_index;

res = intervalquadrature( qdeg, f) * grid.EL( tria_index, face_index);

end