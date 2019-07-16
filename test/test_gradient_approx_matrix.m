function OK = test_gradient_approx_matrix
% function performing a test of the gradient_approx routine by
% generating a simple 2x2 grid on which the gradients over each edge
% can be pre-calculated manually (c.f. grad_correct vector), and
% compared to the output of gradient_approx
% returns 1, if test is OK, 0 otherwise

% This program is open source.  For license terms, see the COPYING file.
%
% --------------------------------------------------------------------
% ATTRIBUTION NOTICE:
% This product includes software developed for the RBmatlab project at
% (C) Universities of Stuttgart and Münster, Germany.
%
% RBmatlab is a MATLAB software package for model reduction with an
% emphasis on Reduced Basis Methods. The project is maintained by
% M. Dihlmann, M. Drohmann, B. Haasdonk, M. Ohlberger and M. Schaefer.
% For Online Documentation and Download we refer to www.morepas.org.
% --------------------------------------------------------------------


% Martin Drohmann 21.05.2008

gradient_approx_matrix_common_settings

U = [1, 0, 1, 0];

grad_correct = [-2; 0;2; 0;-2; 0;2; 0];
OK = check_equal(grad_correct, grad_test(params,model_data,U,[],1),1);
grad_correct = [-1; 0;1; 0;-1; 0;1; 0];
OK = check_equal(grad_correct, grad_test(params,model_data,U,[],2),2) & OK;
grad_correct = [0; 0;-2; 0;0; 0;-2; 0];
OK = check_equal(grad_correct, grad_test(params,model_data,U,[],3),3) & OK;
grad_correct = [-1; 0;1; 0;-1; 0;1; 0];
OK = check_equal(grad_correct, grad_test(params,model_data,U,[],4),4) & OK;
%  grad = grad_test2(U, model_data, params, i)

U = [1, 1, 0, 0];

grad_correct = [0; -1; 0; 0; 0; -1; 2; 0];
OK = check_equal(grad_correct, grad_test(params,model_data,U,[],1),1) & OK;
grad_correct = [-1; -2;1; -2;-2; 0;2; 0];
OK = check_equal(grad_correct, grad_test(params,model_data,U,[],2),2) & OK;
grad_correct = [0; 0;0; -1;-2; 0;0; -1];
OK = check_equal(grad_correct, grad_test(params,model_data,U,[],3),3) & OK;
grad_correct = [0; 0;0; 0;-1; -2;1; -2];
OK = check_equal(grad_correct, grad_test(params,model_data,U,[],4),4) & OK;


f = @(x) sin(x(:,1)) .* cos(x(:,2));
fx = @(x) cos(x(:,1)) .* cos(x(:,2));
fy = @(x) (-1)*sin(x(:,1)) .* sin(x(:,2));

OK = check_convergence(params,f,fx,fy,1) & OK;

end

function [ret] = grad_test(params, model_data,U,NU_ind,edge)
  [grad,bdir] = gradient_approx_matrix(params, model_data, NU_ind, edge);
  ret = grad * U' + bdir;
end

function [OK] = check_equal(correct, testval,edge)
  OK = 1;
  if ~all(all(correct == testval))
    disp(['gradient over edge ', num2str(edge), ' incorrect!']);
    disp(['expected: ', mat2str(correct), ...
          ' got: ',      mat2str(testval)]);
    OK = 0;
  end

end


function [OK] = check_convergence(params,f,fx,fy,edge)
  OK = 1;
  params.xnumintervals = 4;
  params.ynumintervals = 4;
  params.bnd_rect_index = [-2, -2];

  maxerr = zeros(1,4);

  for refstep=1:4
    model_data = nonlin_evol_gen_model_data(params);

    grid = model_data.grid;

    get_enbi(grid);

    dx = grid.DC(1,1);
    dy = grid.DC(1,2);
    ie = find(grid.ECX(:,edge) > dx & grid.ECX(:,edge) < 1-dx...
              & grid.ECY(:,edge) > dy & grid.ECY(:,edge) < 1-dy); % inner edges

    U = f([grid.CX,grid.CY])';
    g1 = reshape([fx([grid.ECX(ie,edge),grid.ECY(ie,edge)])';...
                  fy([grid.ECX(ie,edge),grid.ECY(ie,edge)])'],...
                 2*length(ie), 1);

    recons = grad_test(params,model_data,U,[],edge);
    nie = reshape([2*ie'-1;2*ie'],2*length(ie),1);
    recons = recons(nie,:);

    maxerr(refstep) = max(g1 - recons);

    params.xnumintervals = params.xnumintervals*2;
    params.ynumintervals = params.ynumintervals*2;

  end
  eocrefstep = log(maxerr(1:end-1)./maxerr(2:end))/log(2);
  if min(eocrefstep) < 1.5
    disp(['gradient over edge ', num2str(edge), ' incorrect!']);
    disp(['EOC sequence has too small entries(<1.5): ', mat2str(eocrefstep)]);
    OK = 0;
  end
end


%| \docupdate
%
