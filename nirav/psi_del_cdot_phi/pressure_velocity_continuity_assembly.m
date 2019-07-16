function [ res ] = pressure_velocity_continuity_assembly...
    ( params, paramsP, grid, qdeg)
%PRESSURE_VELOCITY_CONTINUITY_ASSEMBLY Summary of this function goes here:
%L_2 scalar product of (\psi,\nabla \cdot \phi) where \phi is velocity basis function
%and \psi is pressure basis function
%   Detailed explanation goes here:TODO

if nargin == 3
    qdeg = params.qdeg;
end

res = zeros(paramsP.ndofs_per_element * grid.nelements,...
    params.ndofs_per_element * grid.nelements);
gids_pressure = ldg_global_dof_index(paramsP,grid);
gids_velocity = ldg_global_dof_index(params,grid);

for k=1:1:grid.nelements
    A = pressure_velocity_continuity_integral...
        ( params, paramsP, k, grid, qdeg);
    res(gids_pressure(k,:),gids_velocity(k,:))=...
        res(gids_pressure(k,:),gids_velocity(k,:)) + A;
end
if params.show_sparsity == true
    figure()
    spy(res)
    title('\psi_i,\nabla \cdot \phi_j')
    disp('Observe all graphs of \nabla \cdot \phi_j')
    pause()
end
end