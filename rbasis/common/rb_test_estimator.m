function test_err = rb_test_estimator(model, ...
                                      reduced_data,M_test)
%function test_err = rb_test_estimator(reduced_data,M_test)
%
% function determining the test-error-estimators for the given set of
% vectors mu (columns in M_test) of the RB simulation with corresponding 
% RB set. 

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
%rb_simulation = model.rb_simulation;

nmus = size(M_test,2);  

test_err = zeros(nmus,1);  
parfor i = 1:nmus

  tmodel = model;
  if tmodel.verbose>=5
    fprintf('.');
  end;
  tmodel = tmodel.set_mu(model,M_test(:,i),true);
  sim_data = rb_simulation(tmodel,reduced_data);
  test_err(i) = tmodel.get_estimator_from_sim_data(sim_data);

end;
fprintf('\n');

