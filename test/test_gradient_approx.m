function OK = test_gradient_approx
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

  OK = 1;

  gradient_approx_matrix_common_settings

  U = [1, 0, 1, 0];

  grad_correct = [-2 0;2 0;-2 0;2 0];
  if ~all(all(grad_correct == gradient_approx(params,model_data,U,[],1)))
    disp('x gradient over edge 1 incorrect!');
    OK = 0;
  end
  grad_correct = [-1 0;1 0;-1 0;1 0];
  if ~all(all(grad_correct == gradient_approx(params,model_data,U,[],2)))
    disp('x gradient over edge 2 incorrect!');
    OK = 0;
  end
  grad_correct = [0 0;-2 0;0 0;-2 0];
  if ~all(all(grad_correct == gradient_approx(params,model_data,U,[],3)))
    disp('x gradient over edge 3 incorrect!');
    OK = 0;
  end
  grad_correct = [-1 0;1 0;-1 0;1 0];
  if ~all(all(grad_correct == gradient_approx(params,model_data,U,[],4)))
    disp('x gradient over edge 4 incorrect!');
    OK = 0;
  end
  %  grad = gradient_approx2(U, model_data, params, i)

  U = [1, 1, 0, 0];

  grad_correct = [0 -1; 0 0; 0 -1; 2 0];
  if ~all(all(grad_correct == gradient_approx(params,model_data,U,[],1)))
    disp('y gradient over edge 1 incorrect!');
    OK = 0;
  end
  grad_correct = [-1 -2;1 -2;-2 0;2 0];
  if ~all(all(grad_correct == gradient_approx(params,model_data,U,[],2)))
    disp('y gradient over edge 2 incorrect!');
    OK = 0;
  end
  grad_correct = [0 0;0 -1;-2 0;0 -1];
  if ~all(all(grad_correct == gradient_approx(params,model_data,U,[],3)))
    disp('y gradient over edge 3 incorrect!');
    OK = 0;
  end
  grad_correct = [0 0;0 0;-1 -2;1 -2];
  if ~all(all(grad_correct == gradient_approx(params,model_data,U,[],4)))
    disp('y gradient over edge 4 incorrect!');
    OK = 0;
  end

  f = @(x) sin(x(:,1)) .* cos(x(:,2));
  fx = @(x) cos(x(:,1)) .* cos(x(:,2));
  fy = @(x) (-1)*sin(x(:,1)) .* sin(x(:,2));

  for edge = 1:4
    OK = check_convergence(params,f,fx,fy,edge) & OK;
  end
end

function [ret] = grad_test(params, model_data,U,NU_ind,edge)
  ret = gradient_approx(params, model_data, U, NU_ind, edge);
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
    ie = grid.ECX(:,edge) > dx & grid.ECX(:,edge) < 1-dx...
       & grid.ECY(:,edge) > dy & grid.ECY(:,edge) < 1-dy; % inner edges
    U = f([grid.CX,grid.CY])';
    g1 = [fx([grid.ECX(ie,edge),grid.ECY(ie,edge)]),...
          fy([grid.ECX(ie,edge),grid.ECY(ie,edge)])];

    recons = grad_test(params,model_data,U,[],edge);
    recons = recons(ie,:);

    maxerr(refstep) = max(max(g1 - recons));

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
