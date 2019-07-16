function [ res ] = pressure_velocity_continuity_local( lcoord, params, paramsP,...
    k, grid)
%PRESSURE_VELOCITY_CONTINUIUTY_LOCAL Summary of this function goes here:
% L_2 scalar product
%   Detailed explanation goes here: TODO

JIT=[grid.JIT(k,:,1)',grid.JIT(k,:,2)'];

pressure_basis=ldg_evaluate_basis(lcoord,paramsP);

velocity_basis_derivative=ldg_evaluate_basis_derivative(lcoord,params);

res=zeros(paramsP.ndofs_per_element,params.ndofs_per_element);


for i=1:1:paramsP.ndofs_per_element
    for j=1:1:params.ndofs_per_element
        a = rem(j,params.dimrange);
        if a==0
            a=params.dimrange;
        end
        temp = velocity_basis_derivative{j}*JIT';
        res(i,j)=res(i,j)+pressure_basis(i,:)*temp(a,a);
    end
end


% for b=1:1:paramsP.dimrange
%     for i=b:paramsP.dimrange:paramsP.ndofs_per_element
%         for a=1:1:params.dimrange
%             for j=a:params.dimrange:params.ndofs_per_element
%                 %temp = JIT*sum(velocity_basis_derivative{j},1)';
%                 temp = velocity_basis_derivative{j}*JIT';
%                 res(i,j)=res(i,j)+pressure_basis(i,:)*temp(a);
%             end
%         end
%     end
% end
end