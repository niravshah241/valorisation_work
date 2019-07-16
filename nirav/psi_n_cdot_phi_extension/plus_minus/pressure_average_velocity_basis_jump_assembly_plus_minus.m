function [res] = pressure_average_velocity_basis_jump_assembly_plus_minus...
    (grid,params,paramsP)

[ tria_index_internal,local_edge_index_internal,a_index_internal, b_index_internal,...
    local_vertex_index_internal] = tria_edge_index_internal( grid );
[ tria_index_dirichlet,local_edge_index_dirichlet,a_index_dirichlet, b_index_dirichlet,...
    local_vertex_index_dirichlet] = tria_edge_index_dirichlet( grid );

tria_index = [tria_index_internal, tria_index_dirichlet];
local_edge_vertex = [local_edge_index_internal, local_edge_index_dirichlet];
local_vertex_index = [local_vertex_index_internal, local_vertex_index_dirichlet];

gids_velocity = ldg_global_dof_index(params, grid);
gids_pressure = ldg_global_dof_index(paramsP, grid);

% Initialising matrx for each term

res_plus_minus_11 = sparse(zeros(paramsP.ndofs_per_element * ...
    grid.nelements, params.ndofs_per_element * ...
    grid.nelements));
res_plus_minus_12 = sparse(zeros(paramsP.ndofs_per_element * ...
    grid.nelements, params.ndofs_per_element * ...
    grid.nelements));
res_plus_minus_21 = sparse(zeros(paramsP.ndofs_per_element * ...
    grid.nelements, params.ndofs_per_element * ...
    grid.nelements));
res_plus_minus_22 = sparse(zeros(paramsP.ndofs_per_element * ...
    grid.nelements, params.ndofs_per_element * ...
    grid.nelements));

%Assembly of four matrices

for i=1:1:length(tria_index)
    tria_index_neighbour = grid.NBI(tria_index(i),local_vertex_index(i));
    if tria_index_neighbour>0
        res_plus_minus_11(gids_pressure(tria_index(i),:),...
            gids_velocity(tria_index_neighbour,:))... 
            = res_plus_minus_11(gids_pressure(tria_index(i),:),...
            gids_velocity(tria_index_neighbour,:))+...
            pressure_average_velocity_basis_jump_int_plus_int_minus_11...
            (grid,params,paramsP,tria_index(i),local_vertex_index(i));
    end
end

for i=1:1:length(tria_index)
    tria_index_neighbour = grid.NBI(tria_index(i),local_vertex_index(i));
    if tria_index_neighbour>0
        res_plus_minus_12(gids_pressure(tria_index(i),:),...
            gids_velocity(tria_index_neighbour,:))... 
            = res_plus_minus_12(gids_pressure(tria_index(i),:),...
            gids_velocity(tria_index_neighbour,:))+...
            pressure_average_velocity_basis_jump_int_plus_int_minus_12...
            (grid,params,paramsP,tria_index(i),local_vertex_index(i));
    end
end

for i=1:1:length(tria_index)
    tria_index_neighbour = grid.NBI(tria_index(i),local_vertex_index(i));
    if tria_index_neighbour>0
        res_plus_minus_21(gids_pressure(tria_index(i),:),...
            gids_velocity(tria_index_neighbour,:))... 
            = res_plus_minus_21(gids_pressure(tria_index(i),:),...
            gids_velocity(tria_index_neighbour,:))+...
            pressure_average_velocity_basis_jump_int_plus_int_minus_21...
            (grid,params,paramsP,tria_index(i),local_vertex_index(i));
    end
end

for i=1:1:length(tria_index)
    tria_index_neighbour = grid.NBI(tria_index(i),local_vertex_index(i));
    if tria_index_neighbour>0
        res_plus_minus_22(gids_pressure(tria_index(i),:),...
            gids_velocity(tria_index_neighbour,:))... 
            = res_plus_minus_22(gids_pressure(tria_index(i),:),...
            gids_velocity(tria_index_neighbour,:))+...
            pressure_average_velocity_basis_jump_int_plus_int_minus_22...
            (grid,params,paramsP,tria_index(i),local_vertex_index(i));
    end
end

% Plotting and storgage in output
close all
res.plus_minus_11 = res_plus_minus_11;
res.plus_minus_12 = res_plus_minus_12;
res.plus_minus_21 = res_plus_minus_21;
res.plus_minus_22 = res_plus_minus_22;
res.res = res_plus_minus_11 + res_plus_minus_12 + ...
    res_plus_minus_21 + res_plus_minus_22;
if params.show_sparsity == true
    figure()
    spy(res.plus_minus_11)
    title('Spy of ({\psi^+},[n^- \cdot \phi^-])_11')
    figure()
    spy(res.plus_minus_12)
    title('Spy of ({\psi^+},[n^- \cdot \phi^-])_12')
    figure()
    spy(res.plus_minus_21)
    title('Spy of ({\psi^+},[n^- \cdot \phi^-])_21')
    figure()
    spy(res.plus_minus_22)
    title('Spy of ({\psi^+},[n^- \cdot \phi^-])_22')
    figure()
    spy(res.res)
    title('Spy of ({\psi},[n \cdot \phi])_+-')
    sprintf('Observe all graphs %s','{\psi},[n \cdot \phi]')
    pause();
    close all
end

end