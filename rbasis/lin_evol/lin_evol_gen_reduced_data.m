function reduced_data = lin_evol_gen_reduced_data(model,detailed_data)
%function reduced_data = lin_evol_gen_reduced_data(model,detailed_data)
%
% method which produces reduced_data, which is the data, that will be passed to
% an online-algorithm. Therefore, no quantities dependent on the high-dimension
% H may be included here. Neither may online-data include parameter-dependent
% mu-quantities. So no complete grid or detailed solutions or reduced basis
% vectors may be stored here. So online data is produced in the offline stage,
% but may be used in online-stages. So the computation time may depend on H,
% but the results may not depend on this complexity.
%    
% allowed dependency of generated data: Nmax
% not allowed dependency of data: H
% allowed dependency of computation: Nmax, H
% Unknown at this stage: mu, N
%
% Required fields of model:
%   none

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


% Bernard Haasdonk 16.5.2007

model.decomp_mode = 1;
reduced_data.a0 = rb_init_values(model, detailed_data);

% assuming that components do not change in time, so wlg t = 0!!!
%model = set_time(model,0);
[reduced_data.LL_I, reduced_data.LL_E, reduced_data.bb, ...
 reduced_data.K_II, reduced_data.K_IE, reduced_data.K_EE, ... 
 reduced_data.m_I, reduced_data.m_E, reduced_data.m ] = rb_operators(model, detailed_data);
if isfield(model,'name_output_functional')
  reduced_data = rb_output_functional_prep(model,reduced_data,detailed_data);
end;

reduced_data.N = model.get_rb_size(model,detailed_data);



