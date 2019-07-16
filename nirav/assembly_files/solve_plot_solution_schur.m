function [ params, paramsP, achieved_residual_tol_schur,...
    cholesky_error, cholesky_error_inverse] =...
    solve_plot_solution_schur( params, paramsP, grid, rhs, stifness_matrix)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

disp('Entering into solve (Schur)')

L = chol(full(params.bilinear_side),'lower');
inverse_matrix = L'\(L\eye(params.ndofs));

A = params.bilinear_side;

B = params.bilinear_side_pressure_terms;

C = params.bilinear_side_pressure_terms';

a = params.linear_side;

b = params.rhs_continuity;

paramsP.dofs = (-C*inverse_matrix*B)\(b-(C*inverse_matrix*a));
%check formula
params.dofs = (inverse_matrix*a-inverse_matrix*B*paramsP.dofs);
%check formula

dofs = [params.dofs;paramsP.dofs];

achieved_residual_tol_schur = norm(rhs - stifness_matrix*dofs)/norm(rhs);

if nargout>3
    
    cholesky_error = norm(full(params.bilinear_side) - L * L',2);
    
    disp(['Error in cholesky decomposition is ',num2str(cholesky_error)]);
    
    cholesky_error_inverse = norm(inverse_matrix-inv(full(params.bilinear_side)),2);
    
    disp(['Error in inverse with cholesky decomposition is ',...
        num2str(cholesky_error_inverse)]);
    
end

% disp('entering into plotting Degrees of Freedom (Schur)')

% for i=1:1:params.dimrange
%     figure()
%     axis equal
%     [scalar_dofs, scalar_df_info] = ldg_scalar_component(params,i);
%     sdf = ldgdiscfunc(scalar_dofs,scalar_df_info);
%     disp(['Plotting ',num2str(i),' degree of freedom (Schur)'])
%     %subplot(params.dimrange,1,i)
%     %title(['Velocity degree of freedom number ',num2str(i)])
%     if i==1
%         title(['Velocity in x direction (Schur)'])
%     else
%         title(['Velocity in y direction (Schur)'])
%     end
%     axis equal
%     axis tight
%     ldg_plot(sdf,grid,params);
%     plot(grid);
% end
% 
% for i=1:1:paramsP.dimrange
%     figure()
%     [scalar_dofs, scalar_df_info] = ldg_scalar_component(paramsP,i);
%     sdf = ldgdiscfunc(scalar_dofs,scalar_df_info);
%     disp(['Plotting ',num2str(i),' degree of freedom (for pressure)'])
%     %subplot(paramsP.dimrange,1,i)
%     title('Pressure (Schur)')
%     axis equal
%     axis tight
%     ldg_plot(sdf,grid,paramsP);
%     plot(grid);
% end

end