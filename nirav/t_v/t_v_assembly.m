function [res ] = t_v_assembly( params, paramsP, grid, qdeg )
%T_V_ASSEMBLY Summary of this function goes here
%   Detailed explanation goes here

[tria_index,local_edge_index,a_index, b_index, local_vertex_index] =...
    tria_edge_index_neumann( grid );

res=zeros(params.ndofs,1);
gids = ldg_global_dof_index(params,grid);


for i=1:1:length(tria_index)
    res(gids(tria_index(i),:),1) = res(gids(tria_index(i),:),1) + ...
        t_v_integral( params, paramsP, grid, tria_index(i), ...
        local_vertex_index(i), qdeg);
end

end