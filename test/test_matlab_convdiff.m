function OK = test_matlab_convdiff
% function performing a test of the elementary convdiff_model
% routines. In particular simulation without diffusion must be
% 10 times faster!!
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


% Bernard Haasdonk 21.7.2009

OK = 1;

model      = convdiff_model([]);
model_data = gen_model_data(model);
% without diffusion:
model = model.set_mu(model, [1,0,0]);
tic;
detailed_simulation(model,model_data);
t1 = toc;
% with diffusion
model = model.set_mu(model, [1, 0, 5e-8]);
tic;
detailed_simulation(model,model_data);
t2 = toc;

if t2/t1<5
  error('no! diffusion computation should be 10-20 times faster !!')
  OK = 0;
end;


%| \docupdate 
