function [max_test_errs, max_mu_index, ...
	  min_test_errs, min_mu_index] = rb_test_convergence(detailed_data,model)
%function [max_test_errs, max_mu_index,min_test_errs, min_mu_index] =
%                              rb_test_convergence(detailed_data,model)
%
% function determining the maximum and minimum test-error (linfty-l2 or
% estimator) for the given set of vectors mu (columns in 
% detailed_data.RB_info.M_test) of the RB simulation with corresponding
% RB set. A convergence test is
% performed by performing this maximum and minimum
% detection for all numbers of RB from 1 to the complete set.
% the output is max_test_errs(i): the maximum test-quantity on the given
% mu-subset for reduced basis RB(:,1:i). The number of the mu, which
% incurs this maximum error is returned in
% max_mu_index(i). Similarly for min_test_errs and min_mu_index.
%
% Required fields of model:
%        RB_error_indicator: 'error' or 'estimator'
%        RB_detailed_test_savepath : in case of 'error' this path
%               either contains the test-samples or they are
%               generated.
%        error_algorithm : algorithm computing the true error in
%               case of 'error' mode
%        test_N_samples : (optional) number of N samples, which are tested,
%                   i.e. value 11 for a basis of size N=21 will give 11
%                   test-results for the values N = 1,3,5,7,9,11,13,15,17,19,21
%                   an equidistant sampling of the N-interval is
%                   realized. If not specified, all numbers 1:N are tested

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

% Bernard Haasdonk 2.8.2006
  
nmus = size(detailed_data.RB_info.M_test,2);  

nRBs = size(detailed_data.RB,2);

max_test_errs = NaN(nRBs,1);  
max_mu_index = NaN(nRBs,1);  
min_test_errs = NaN(nRBs,1);  
min_mu_index = NaN(nRBs,1);  

% compute offline data once for complete RB set
model.Nmax = nRBs;
offline_data  = rb_offline_prep(detailed_data,model);

if isequal(model.RB_error_indicator,'error')
  savepath = model.RB_detailed_test_savepath;
  
  if isempty(savepath)
    error('savepath must be provided in case of error as target value!');
  end;
  save_detailed_simulations(detailed_data.RB_info.M_test,...
			    savepath,...
			    detailed_data.grid,...
			    model);
end;

if ~isfield(model,'test_N_samples') || ...
      isempty(model.test_N_samples) % default: all numbers from N to 1
  Nsamples = (nRBs:-1:1);
else % compute equidistributed sample-number
  Nsamples = ceil((0:(model.test_N_samples-1))*...
		  (nRBs-1)/(model.test_N_samples-1)+ 0.5);
  Nsamples = Nsamples(end:-1:1);
%  keyboard;
end;

for nRB = Nsamples
  disp(['testing for nrb = ',num2str(nRB)]);
  model.N = nRB;
  reduced_data = rb_online_prep(offline_data,model);
  
  for i = 1:nmus
    model = model.set_mu(detailed_data.RB_info.M_test(:,i),model);
    simulation_data = rb_simulation(reduced_data,model);
    
    if isequal(model.RB_error_indicator,'error')
      sim_data = load_detailed_simulation(i, ...
                                          savepath,...
                                          model);
      U = sim_data.U;
      Uappr = rb_reconstruction(detailed_data,simulation_data);
      errs = feval(model.error_algorithm,U,Uappr,detailed_data.grid,model);
      err = errs(end);
      
    elseif isequal(model.RB_error_indicator, 'estimator')
      err = simulation_data.Delta(end);
    else
      error('error indicator type unknown');
    end;
    test_err(i) = err;   
  end;
  
  [max_test_errs(nRB), max_mu_index(nRB)] = max(test_err);
  [min_test_errs(nRB), min_mu_index(nRB)] = min(test_err);
end;
