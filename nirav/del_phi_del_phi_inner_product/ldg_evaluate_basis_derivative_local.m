function res = ldg_evaluate_basis_derivative_local...
    ( lcoord, params, k, grid)
%LDG_EVALUATE_BASIS_DERIVATIVE_LOCAL: assembly of \phi:\phi where \phi is
%global basis function for velocity. \nabla \phi = JIT * \hat{\nabla \phi} where
%\hat{\phi} is local basis function for velocity
%Params needs fields: dimrange, pdeg, ndofs, ndofs_per_element
%k is element number and grid should have field JIT (Jacobian Inverse Transposed)

basis_derivative = ldg_evaluate_basis_derivative(lcoord,params);

JIT = [grid.JIT(k,:,1)',grid.JIT(k,:,2)']; % Jacobian Inverse Transpose

%Crate global derivative from local derivative and make matrix of global
%derivatives (to avoid for loop)

local_derivative_matrix = cell2mat(basis_derivative);

local_derivative_matrix_reshape = reshape(local_derivative_matrix, ...
    params.dimrange, 2, params.ndofs_per_element);%2 is for 2-D grid

local_derivative_matrix_JIT = reshape(local_derivative_matrix_reshape,...
    params.dimrange * 2, 1, params.ndofs_per_element);

local_derivative_matrix_JIT = reshape(local_derivative_matrix_JIT,...
    params.dimrange * 2, params.ndofs_per_element, 1);

res.local_derivative_matrix_JIT = local_derivative_matrix_JIT;

JIT_modified = zeros(4);

JIT_modified(1,1) = JIT(1,1);
JIT_modified(3,1) = JIT(1,2);
JIT_modified(2,2) = JIT(1,1);
JIT_modified(4,2) = JIT(1,2);
JIT_modified(1,3) = JIT(2,1);
JIT_modified(3,3) = JIT(2,2);
JIT_modified(2,4) = JIT(2,1);
JIT_modified(4,4) = JIT(2,2);

global_derivative_matrix_reshaped = local_derivative_matrix_JIT' * ...
    JIT_modified;

% global_derivative_matrix = reshape(global_derivative_matrix_reshaped,...
%     params.dimrange,2,params.ndofs_per_element);
% 
% global_derivative_matrix = reshape(global_derivative_matrix,2*params.dimrange,...
%     1,params.ndofs_per_element);
% 
% global_derivative_matrix = reshape(global_derivative_matrix,2*params.dimrange,...
%     params.ndofs_per_element,1)';

% Now in global derivative matrix each row means one \nabla \phi

res = global_derivative_matrix_reshaped * ...
    global_derivative_matrix_reshaped'; 

res = 1/2 * ( res + res');

%res=zeros(params.ndofs_per_element);

% for a=1:1:params.dimrange
%     for i=a:params.dimrange:params.ndofs_per_element
%         for j=a:params.dimrange:params.ndofs_per_element
%             res(i,j)=sum(basis_derivative{i},1)*JIT'*JIT*...
%                 sum(basis_derivative{j},1)';...
%         end
%     end
% end

%res = (res+res')/2;

end