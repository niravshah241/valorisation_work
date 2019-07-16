function [ res ] = pressure_velocity_continuity_assembly_extension...
    ( params, paramsP, grid, qdeg)
%PRESSURE_VELOCITY_CONTINUITY_ASSEMBLY Summary of this function goes here:
%L_2 scalar product of (\psi,\nabla \cdot \phi) where \phi is velocity basis function
%and \psi is pressure basis function
%   Detailed explanation goes here:TODO

if nargin == 3
    qdeg = params.qdeg;
end

res_11 = sparse(zeros(paramsP.ndofs_per_element*grid.nelements,...
    params.ndofs_per_element*grid.nelements));
res_12 = sparse(zeros(paramsP.ndofs_per_element*grid.nelements,...
    params.ndofs_per_element*grid.nelements));
res_21 = sparse(zeros(paramsP.ndofs_per_element*grid.nelements,...
    params.ndofs_per_element*grid.nelements));
res_22 = sparse(zeros(paramsP.ndofs_per_element*grid.nelements,...
    params.ndofs_per_element*grid.nelements));

gids_pressure = ldg_global_dof_index(paramsP,grid);
gids_velocity = ldg_global_dof_index(params,grid);

for k=1:1:grid.nelements
    A_11 = pressure_velocity_continuity_integral_11_extension...
        ( params, paramsP, k, grid);
    res_11(gids_pressure(k,:),gids_velocity(k,:))=...
        res_11(gids_pressure(k,:),gids_velocity(k,:)) + A_11;
    A_12 = pressure_velocity_continuity_integral_12_extension...
        ( params, paramsP, k, grid);
    res_12(gids_pressure(k,:),gids_velocity(k,:))=...
        res_12(gids_pressure(k,:),gids_velocity(k,:)) + A_12;
    A_21 = pressure_velocity_continuity_integral_21_extension...
        ( params, paramsP, k, grid);
    res_21(gids_pressure(k,:),gids_velocity(k,:))=...
        res_21(gids_pressure(k,:),gids_velocity(k,:)) + A_21;
    A_22 = pressure_velocity_continuity_integral_22_extension...
        ( params, paramsP, k, grid);
    res_22(gids_pressure(k,:),gids_velocity(k,:))=...
        res_22(gids_pressure(k,:),gids_velocity(k,:)) + A_22;
end

res{1} = res_11;
res{2} = res_12;
res{3} = res_21;
res{4} = res_22;

end