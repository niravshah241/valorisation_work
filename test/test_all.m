function OK = test_all
% function OK = test_all
%
% function calling all tests and printing diagnostics.  
% should be called after major changes. 
%
% the called tests should output some information, that can be
% interpreted. The called tests should not break with an error, such that
% this routine will consecutively perform all tests, even if a failed
% test occurs.

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
    
test_names = {'test_matlab_convdiff','test_fail','test_gradient_approx',...
              'test_gradient_approx_matrix', ...
              'test_lebesgue', ...
              'test_rb_local_ext', ...
              'test_orthonormalize',...
	      'test_rb_lin_evol', 'test_rb_basisgen', ...
              'test_power_vector2',...
	      'test_ldgfunc','test_ldg_derivative','test_ldg_orthogonality'};

OKs = zeros(1,length(test_names));

for i = 1:length(test_names);
  disp('-----------------------------------------------------------')
  disp(['testing ', test_names{i}]);
  try
    OKs(i) = feval(test_names{i});
  catch
    warning(['test ',test_names{i},' failed!']);
    OKs(i) = 0;
  end;
  if OKs(i)==0
    disp(' => failed')
  else    
    disp(' => success')
  end;
end;

disp('-----------------------------------------------------------')
disp('test results:');
resstr = {'failed', 'OK'};
for i = 1:length(test_names);
  disp([test_names{i},' => ' resstr{OKs(i)+1}]);
end;
disp('-----------------------------------------------------------')

OK = all(OKs);



%| \docupdate 
