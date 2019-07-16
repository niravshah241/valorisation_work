function [ res ] = q_n_u_d_assembly( params, paramsP, grid, qdeg )
%Q_N_U_D_ASSEMBLY Summary of this function goes here
%   Detailed explanation goes here

[ tria_index,local_edge_index,a_index, b_index, local_vertex_index]...
    = tria_edge_index_dirichlet( grid );

res=zeros(paramsP.ndofs,1);

gids_pressure = ldg_global_dof_index(paramsP, grid);

for i=1:1:length(tria_index)
    res(gids_pressure(tria_index(i),:))=res(gids_pressure(tria_index(i),:)) +... 
        q_n_u_d_integral( params, paramsP,grid,tria_index(i),local_vertex_index(i), qdeg);
end

end