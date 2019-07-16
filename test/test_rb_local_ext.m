function OK = test_rb_local_ext
% function performing a test wether all necessary cells surrounding an
% interpolation point are selected during the reduced_data generation phase.
%
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

gradient_approx_matrix_common_settings;

%detailed_data = structcpy([], model_data);

params.xnumintervals = 5;
params.ynumintervals = 5;

params.stencil_mode = 'vertex';

detailed_data = nonlin_evol_gen_model_data(params);

% we define, there is only one interpolation point at index 8
detailed_data.TM{1} = 8;
% create some dummy reduced basis / collateral basis
detailed_data.QM{1} = ones(25,1);
detailed_data.RB    = (1:25)';
detailed_data.BM{1} = 1;

detailed_data.implicit_crb_index = 1;
detailed_data.explicit_crb_index = 1;

% the points we expect to be selected
ref_RB = [2:4,7:9,12:14]';
reduced_data = nonlin_evol_gen_reduced_data(params, detailed_data);

if length(reduced_data.RB_local_ext{1}) ~= length(ref_RB) ...
    || ~all(all(reduced_data.RB_local_ext{1} == ref_RB))
  disp('local gridpart construction is faulty!');
  OK = 0;
end

end

%| \docupdate 
