function OK = test_orthonormalize
% function performing a test of the orthonormalization routine by
% generating a time-sequence of data and orthogonalizing this.
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


% Bernard Haasdonk

OK = 1;

load('test_rb_lin_evol_data');
model_data = gen_model_data(model);
sim_data = detailed_simulation(model,model_data);
A = model_data.W;
UON = orthonormalize(sim_data.U,A);
B = UON' * A * UON;

E = eye(size(B));
i_zero = find(E==0);
i_nonzero = find(E~=0);

if (find(B>0.9999)~=i_nonzero)
  disp('gram-matrix not the identity!');
  OK = 0;
end;

if (find(B<0.0001)~=i_zero)
  disp('gram-matrix not the identity!');
  OK = 0;
end;

%pcolor(B);

%| \docupdate 
