function test_err = rb_test_indicator(model,...
                                      detailed_data,reduced_data,...
                                      M_test,savepath)
%function test_err = rb_test_indicator(model,[detailed_data],[reduced_data]
%                                      M_test,[savepath])
%
% function determining the test-error-indicators for the given set of
% vectors mu (columns in M_test) of the RB simulation with corresponding
% RB set. Either the true errors or the error-estimators are
% determined. If the true errors are wanted, these are produced
% into or read from the directory savepath. The offline_data is
% generated, if empty. The detailed data can be ommited if
% offline_data is available, and no true 'error' is wanted
%
% Required fields of model:
%    RB_error_indicator: 'error' or 'estimator'
%    error_algorithm : in case of 'error' mode, this algorithm
%           is called for computing the error between a U and Uappr

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

% Bernard Haasdonk 27.5.2007

if isempty(reduced_data)
  reduced_data = gen_reduced_data(model,detailed_data);
end;

if ~isfield(model,'N')
    if isfield(detailed_data,'N')
        model.N=detailed_data.N;
    else
        model.N = size(model.get_rb_from_detailed_data(detailed_data),2);
    end
end;
if isfield(model,'rb_problem_type')
if isequal(model.rb_problem_type,'nonlin_evol') ...
     || isequal(model.rb_problem_type, 'richards')
  if ~isfield(model,'M') || (model.M==-1)
    model.M = max(cellfun(@(X)(size(X,2)), detailed_data.QM));
  end;
end;
end;


reduced_data = model.reduced_data_subset(model, reduced_data);


switch model.RB_error_indicator
 case 'error'

  test_err = rb_test_error(model,...
                           detailed_data, reduced_data, M_test, ...
                           savepath);

 case 'estimator'
  test_err = rb_test_estimator(model,...
                               reduced_data,M_test);

 case 'ei_estimator_test'
  tic;
   test_err.estimator = rb_test_estimator(model, ...
                                          reduced_data,M_test);
  t1 = toc;
  tic;
   test_err.error     = rb_test_error(model, ...
                                      detailed_data, reduced_data, M_test, ...
                                      savepath);
  t2 = toc;
  lambda = test_err.estimator ./ test_err.error;
  disp(['Evaluation of the reliability lambda=Delta_est/Delta_err of a ',...
        'posteriori estimator (min, max, mean):']);
  disp(['min(lambda)=',num2str(min(lambda)), ...
        '  max(lambda)=',num2str(max(lambda)), ...
        '  mean(lambda)=', num2str(mean(lambda))]);
  disp('Evaluation of the runtime gain of posteriori estimator');
  disp(['t_post=',num2str(t1),'   t_err=',num2str(t2)]);
  if(any(test_err.estimator < test_err.error))
    disp('Error: a posteriori estimator underestimates the real error');
    keyboard;
  end
 otherwise
  error('error indicator type unknown');
end;

