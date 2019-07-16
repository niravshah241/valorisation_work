function [res] = n_u_tensor_jump_n_phi_tensor_jump_plus_minus_assembly...
    (grid,params)

[ tria_index_internal,local_edge_index_internal,a_index_internal, ...
    b_index_internal, local_vertex_index_internal] = ...
    tria_edge_index_internal( grid );
[ tria_index_dirichlet,local_edge_index_dirichlet,a_index_dirichlet, ...
    b_index_dirichlet, local_vertex_index_dirichlet] = ...
    tria_edge_index_dirichlet( grid );

tria_index = [tria_index_internal, tria_index_dirichlet];
local_edge_vertex = [local_edge_index_internal, ...
    local_edge_index_dirichlet];
local_vertex_index = [local_vertex_index_internal, ...
    local_vertex_index_dirichlet];

gids_velocity = ldg_global_dof_index(params, grid);

% Initialising matrix for each term

res_plus_minus_11 = sparse(zeros(params.ndofs));
res_plus_minus_12 = sparse(zeros(params.ndofs));
res_plus_minus_21 = sparse(zeros(params.ndofs));
res_plus_minus_22 = sparse(zeros(params.ndofs));

%Assembly of four matrices

for i=1:1:length(tria_index)
    tria_index_neighbour = grid.NBI(tria_index(i),local_vertex_index(i));
    if tria_index_neighbour > 0
        res_plus_minus_11(gids_velocity(tria_index_neighbour,:),...
            gids_velocity(tria_index(i),:)) ...
            = res_plus_minus_11(gids_velocity(tria_index_neighbour,:),...
            gids_velocity(tria_index(i),:)) + ...
            n_u_tensor_jump_n_phi_tensor_jump_plus_minus_integral_11...
            (params,grid,tria_index(i),local_vertex_index(i));
    end
end

for i=1:1:length(tria_index)
    tria_index_neighbour = grid.NBI(tria_index(i),local_vertex_index(i));
    if tria_index_neighbour > 0
        res_plus_minus_12(gids_velocity(tria_index_neighbour,:),...
            gids_velocity(tria_index(i),:)) ...
            = res_plus_minus_12(gids_velocity(tria_index_neighbour,:),...
            gids_velocity(tria_index(i),:)) + ...
            n_u_tensor_jump_n_phi_tensor_jump_plus_minus_integral_12...
            (params,grid,tria_index(i),local_vertex_index(i));
    end
end

for i=1:1:length(tria_index)
    tria_index_neighbour = grid.NBI(tria_index(i),local_vertex_index(i));
    if tria_index_neighbour > 0
        res_plus_minus_21(gids_velocity(tria_index_neighbour,:),...
            gids_velocity(tria_index(i),:)) ...
            = res_plus_minus_21(gids_velocity(tria_index_neighbour,:),...
            gids_velocity(tria_index(i),:)) + ...
            n_u_tensor_jump_n_phi_tensor_jump_plus_minus_integral_21...
            (params,grid,tria_index(i),local_vertex_index(i));
    end
end

for i=1:1:length(tria_index)
    tria_index_neighbour = grid.NBI(tria_index(i),local_vertex_index(i));
    if tria_index_neighbour > 0
        res_plus_minus_22(gids_velocity(tria_index_neighbour,:),...
            gids_velocity(tria_index(i),:)) ...
            = res_plus_minus_22(gids_velocity(tria_index_neighbour,:),...
            gids_velocity(tria_index(i),:)) + ...
            n_u_tensor_jump_n_phi_tensor_jump_plus_minus_integral_22...
            (params,grid,tria_index(i),local_vertex_index(i));
    end
end

% Plotting and storgage in output
%disp('check boundary nodes')
%pause();
res.plus_minus_11 = res_plus_minus_11;
res.plus_minus_12 = res_plus_minus_12;
res.plus_minus_21 = res_plus_minus_21;
res.plus_minus_22 = res_plus_minus_22;
res.res = res_plus_minus_11 + res_plus_minus_12 + res_plus_minus_21 + ...
    res_plus_minus_22;

if params.show_sparsity == true
    figure()
    spy(res.plus_minus_11)
    title('Spy of (([n \cdot \phi]^+ , [n \cdot \phi]^-)_11')
    figure()
    spy(res.plus_minus_12)
    title('Spy of (([n \cdot \phi]^+ , [n \cdot \phi]^-)_12')
    figure()
    spy(res.plus_minus_21)
    title('Spy of (([n \cdot \phi]^+ , [n \cdot \phi]^-)_21')    
    figure()
    spy(res.plus_minus_22)
    title('Spy of (([n \cdot \phi]^+ , [n \cdot \phi]^-)_22')    
    figure()
    spy(res.res)
    title('Spy of ([n \otimes \phi]^+,[n \otimes \phi]^-)')
    sprintf('Observe all graphs %s','[[n \otimes \phi]],[[n \otimes \phi]]')
    pause();
    close all
end

end