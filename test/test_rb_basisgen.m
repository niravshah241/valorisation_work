function OK = test_rb_basisgen
%function OK = test_rb_basisgen
%
% function generating the same basis as stored in
% demo_lin_evol_detailed_data.mat by performing a greedy basis extension
% on a 5x5x5 parameter grid.
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
model_data = gen_model_data(model);

model.detailedfname = fullfile(rbmatlabtemp, 'test_rb_basisgen');

% only generate 20 basis functions, otherwise too long.
params.RB_stop_Nmax = 20;
detailed_data = gen_detailed_data(model, model_data);

if size(detailed_data.RB,2) > 20
  disp('RB basis is too large!')
  OK = 0;
end;

if detailed_data.RB_info.max_err_sequence(20) > 1.2370e-04 * 1.1;
  disp('error(20) is too large !')
  OK = 0;
end;

if detailed_data.RB_info.max_err_sequence(20) < 1.2370e-04 * 0.1;
  disp('error(20) is too small !')
  OK = 0;
end;


%| \docupdate 
