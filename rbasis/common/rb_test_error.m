function test_err = rb_test_error(model,detailed_data,reduced_data,M_test,...
                                  savepath)
%function test_err = rb_test_error(model,detailed_data,reduced_data,M_test,...
%                                 savepath)
%
% function determining the test-errors, i.e. energy norm or l2-error given ty
% 'model.error_algorithm' at end time for the given set of vectors `\mu`
% (columns in 'M_test') of the RB simulation with corresponding RB set.
%
% parameters:
%   M_test: matrix with column vectors of parameter tuples `\mu` for which the
%           error between reduced and detailed simulation shall be computed.
%   savepath: directory path where the detailed simulations shall be saved to
%             with save_detailed_simulations().
%
% required fields of model:
%   error_algorithm: function handle to function that computes the error
%                    `\|u_h - u_{red}\|` in some norm.
%
% optional fields of model:
%   relative error: boolean indicating wether relative or absolute errors shall
%                   be computed. (Default = absolut errors)
%
% return values:
%   test_err: vector of errors at end time

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

% Bernard Haasdonk 6.6.2007

% import model specific methods

nmus = size(M_test,2);
test_err = zeros(nmus,1);

if isempty(savepath)
  error('savepath must be provided in case of error as target value!');
end;
save_detailed_simulations(model, detailed_data, M_test, savepath);

parfor i = 1:nmus
  fprintf('.');
  tmodel = model;

  tmodel = tmodel.set_mu(model, M_test(:,i));
  try
    rb_sim_data = rb_simulation(tmodel,reduced_data);
    sim_data = load_detailed_simulation(i, ...
                                      savepath,...
                                      tmodel);
    rb_sim_data = rb_reconstruction(tmodel,detailed_data,rb_sim_data);

    if isfield(detailed_data,'grid')
      err_par = detailed_data.grid;
    else
      err_par = [];
    end;

    errs = tmodel.error_algorithm(tmodel.get_dofs_from_sim_data(sim_data),...
                                 tmodel.get_dofs_from_sim_data(rb_sim_data),...
                                 err_par, model);
    if isfield(tmodel,'relative_error') && tmodel.relative_error
      errs = errs./tmodel.error_algorithm(tmodel.get_dofs_from_sim_data(sim_data), ...
          0, err_par, tmodel);
    end
  catch e
    errs = inf;
  end
  err = max(errs);
%  err = errs(end);
%  disp('warning, check monotonicity in errs!!!');
  test_err(i) = err;
%  if (i==2) & (mod(size(detailed_data.RB,2),30) == 0)
%    save(['errs',num2str(size(detailed_data.RB,2))],'errs');
%  end;
end;

disp(test_err');

