function [res] = n_u_tensor_jump_n_phi_tensor_jump_assembly...
    (grid,params,qdeg)

if nargin == 2
    qdeg = params.qdeg;
end

[ tria_index_internal,local_edge_index_internal,a_index_internal, b_index_internal,...
    local_vertex_index_internal] = tria_edge_index_internal( grid );
[ tria_index_dirichlet,local_edge_index_dirichlet,a_index_dirichlet, b_index_dirichlet,...
    local_vertex_index_dirichlet] = tria_edge_index_dirichlet( grid );

tria_index = [tria_index_internal, tria_index_dirichlet];
local_edge_vertex = [local_edge_index_internal, local_edge_index_dirichlet];
local_vertex_index = [local_vertex_index_internal, local_vertex_index_dirichlet];

gids_velocity = ldg_global_dof_index(params, grid);

% Initialising matrix for each term

res_plus_plus=zeros(params.ndofs_per_element*grid.nelements);
res_plus_minus=zeros(params.ndofs_per_element*grid.nelements);
res_minus_plus=zeros(params.ndofs_per_element*grid.nelements);
res_minus_minus=zeros(params.ndofs_per_element*grid.nelements);

%Assembly of four matrices

for i=1:1:length(tria_index)
    res_plus_plus(gids_velocity(tria_index(i),:),gids_velocity(tria_index(i),:))...
        =res_plus_plus(gids_velocity(tria_index(i),:),gids_velocity(tria_index(i),:))+...
        n_u_tensor_jump_n_phi_tensor_jump_plus_plus_integral...
        (params,grid,tria_index(i),local_vertex_index(i),qdeg);
end

for i=1:1:length(tria_index)
    tria_index_neighbour = grid.NBI(tria_index(i),local_vertex_index(i));
    if tria_index_neighbour>0
        res_plus_minus(gids_velocity(tria_index_neighbour,:),gids_velocity(tria_index(i),:))...
            =res_plus_minus(gids_velocity(tria_index_neighbour,:),gids_velocity(tria_index(i),:))+...
            n_u_tensor_jump_n_phi_tensor_jump_plus_minus_integral...
            (params,grid,tria_index(i),local_vertex_index(i),qdeg);
    end
end

for i=1:1:length(tria_index)
    tria_index_neighbour = grid.NBI(tria_index(i),local_vertex_index(i));
    if tria_index_neighbour>0
        res_minus_plus(gids_velocity(tria_index(i),:),gids_velocity(tria_index_neighbour,:))...
            =res_minus_plus(gids_velocity(tria_index(i),:),gids_velocity(tria_index_neighbour,:))+...
            n_u_tensor_jump_n_phi_tensor_jump_minus_plus_integral...
            (params,grid,tria_index(i),local_vertex_index(i),qdeg);
    end
end

for i=1:1:length(tria_index)
    tria_index_neighbour = grid.NBI(tria_index(i),local_vertex_index(i));
    if tria_index_neighbour>0
        res_minus_minus(gids_velocity(tria_index_neighbour,:),gids_velocity(tria_index_neighbour,:))...
            =res_minus_minus(gids_velocity(tria_index_neighbour,:),gids_velocity(tria_index_neighbour,:))+...
            n_u_tensor_jump_n_phi_tensor_jump_minus_minus_integral...
            (params,grid,tria_index(i),local_vertex_index(i),qdeg);
    end
end

% Plotting and storgage in output
%disp('check boundary nodes')
%pause();
res.plus_plus = res_plus_plus;
res.plus_minus = res_plus_minus;
res.minus_minus = res_minus_minus;
res.minus_plus = res_minus_plus;
res.res = res_plus_plus + res_plus_minus + res_minus_plus + res_minus_minus;
%res.res = (res.res'+res.res)/2;
if params.show_sparsity == true
    figure()
    spy(res.plus_plus)
    title('Spy of (([n \cdot \phi]^+ , [n \cdot \phi]^+)')
    figure()
    spy(res.plus_minus)
    title('Spy of ([n \cdot \phi]^+ , [n \cdot \phi]^-)')
    figure()
    spy(res.minus_plus)
    title('Spy of ([n \cdot \phi]^- , [n \cdot \phi]^+)')
    figure()
    spy(res.minus_minus)
    title('Spy of ([n \cdot \phi]^- , [n \cdot \phi]^-)')
    figure()
    spy(res.res)
    title('Spy of ([n \otimes \phi],[n \otimes \phi])')
    sprintf('Observe all graphs %s','[[n \otimes \phi]],[[n \otimes \phi]]')
    pause();
    close all
end
%1. titles of graphs
end