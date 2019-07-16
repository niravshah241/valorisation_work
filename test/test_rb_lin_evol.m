function OK = test_rb_lin_evol
%function OK = test_rb_lin_evol
%
% function performing the whole sequence of reduced basis approximation
% and testing, whether the error prediction of the reduced model is
% indeed smaller than the real error.
%
% OK == 1 if test is OK, otherwise 0

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


% Bernard Haasdonk 20.8.2007

  OK = 1;

model = init_model('convdiff_model');
disp('linear evolution problem loaded');

%model_data = gen_model_data(model);

%detailed_data = gen_detailed_data
%save('test_rb_lin_evol_data', 'detailed_data', 'model');
% do not compute the detailed data, but simply load
load('test_rb_lin_evol_data','detailed_data');
detailed_data.grid = construct_grid(model);

model.verbose = 0;
%model.mu=[0,0,0];
%model.RB_stop_Nmax = 100;
%detailed_data = gen_detailed_data(model, model_data);
%save('/tmp/temp.mat','detailed_data','model','model_data');

reduced_data = gen_reduced_data(model, detailed_data);
disp('offline preparation terminated');

model.N = size(detailed_data.RB,2);

reduced_data = extract_reduced_data_subset(model, reduced_data);

disp('online preparation terminated');

mus = {[0,0,0],[0.1,0.1,0.1e-8]};
elim = [5.7e-9, 2.4e-07]; % error limits

for m = 1:length(mus)
  model = model.set_mu(model, mus{m});

  disp('mu is set');

  rb_simulation_data = rb_simulation(model,reduced_data);

  disp('rb simulation terminated');

  rb_simulation_data = rb_reconstruction(model,detailed_data,rb_simulation_data);

  disp('reconstruction finished');

  simulation_data = detailed_simulation(model, detailed_data);

  disp('detailed simulation finished');

  err = model.error_algorithm(simulation_data.U,rb_simulation_data.U,detailed_data.grid,model);

  disp('error computation finished');

  if err(end) > elim(m);
    disp('reduced error is too large!!!');
    OK = 0;
  end;

  if err(end) < elim(m)*0.1;
    disp('reduced error is too small!!!');
    OK = 0;
  end;

  if rb_simulation_data.Delta(end) < err(end)
    disp('rb_prediction smaller than true error!!!');
    OK = 0;
  end;

end;






%| \docupdate 
