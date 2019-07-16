function [ res ] = u_d_v_assembly( grid, params, qdeg)
%U_D_V_ASSEMBLY Summary of this function goes here
%   Detailed explanation goes here

res=zeros(params.ndofs,1);

[ tria_index,local_edge_index,a_index, b_index, local_vertex_index]...
    = tria_edge_index_dirichlet( grid );

res=zeros(params.ndofs_per_element*grid.nelements,1);
gids=ldg_global_dof_index(params,grid);

for i=1:1:length(tria_index)
    res(gids(tria_index(i),:),1)=res(gids(tria_index(i),:),1)+...
        u_d_v_integral(grid, params,tria_index(i), local_vertex_index(i), qdeg);
end

end

