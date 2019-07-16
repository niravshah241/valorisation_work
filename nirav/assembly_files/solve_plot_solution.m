function [ params, paramsP, flag, achieved_residual_tol, actual_iter] = solve_plot_solution...
    ( params, paramsP, grid, rhs, stifness_matrix, required_residual_tol, max_iter)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

disp('Entering into solve')
%  dofs = stifness_matrix\rhs;
%  flag = 'none';
%  achieved_residual_tol = stifness_matrix*dofs - rhs;
[dofs, flag, achieved_residual_tol, actual_iter] = minres(stifness_matrix,rhs,...
   required_residual_tol,max_iter);
% dofs = inv(stifness_matrix)*rhs; %too expensive
%minres,bicgstab(stifness_matrix,rhs,required_residual_tol,max_iter);

params.dofs = dofs(1:params.ndofs);
paramsP.dofs = dofs(params.ndofs+1:params.ndofs+paramsP.ndofs);

disp('entering into plotting Degrees of Freedom')

for i=1:1:params.dimrange
    figure()
    axis equal
    [scalar_dofs, scalar_df_info] = ldg_scalar_component(params,i);
    sdf = ldgdiscfunc(scalar_dofs,scalar_df_info);
    disp(['Plotting ',num2str(i),' degree of freedom'])
    %subplot(params.dimrange,1,i)
    %title(['Velocity degree of freedom number ',num2str(i)])
    if i==1
        title(['Velocity in x direction'])
    else
        title(['Velocity in y direction'])
    end
    axis equal
    axis tight
    ldg_plot(sdf,grid,params);
    plot(grid);
end

for i=1:1:paramsP.dimrange
    figure()
    [scalar_dofs, scalar_df_info] = ldg_scalar_component(paramsP,i);
    sdf = ldgdiscfunc(scalar_dofs,scalar_df_info);
    disp(['Plotting ',num2str(i),' degree of freedom (for pressure)'])
    %subplot(paramsP.dimrange,1,i)
    title('Pressure')
    axis equal
    axis tight
    ldg_plot(sdf,grid,paramsP);
    plot(grid);
end
%pause();
%close all

end