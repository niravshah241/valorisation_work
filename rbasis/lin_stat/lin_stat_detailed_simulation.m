function sim_data = lin_stat_detailed_simulation(model,model_data)
%function sim_data = lin_stat_detailed_simulation(model,model_data)
%
% function performing the detailed simulation of a lin-stat model,
% i.e. the matrix and rhs assembly and solving of the system,
% possibly also output computation

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


% B. Haasdonk 22.2.2011

old_mode = model.decomp_mode;
model.decomp_mode = 0; % complete

sim_data =[];

[A,r] = model.operators(model,model_data);

% solution variable:
uh = femdiscfunc([],model_data.df_info);
uh.dofs = A\r;    

% return results:
sim_data.uh = uh;

% compute output
if model.compute_output_functional
  % the following can be used for any, also nonlinear functionals:
  %sim_data.s =
  %model.output_functional(model,model_data,sim_data.uh);
  % for linear operators, get vector:
  v = model.operators_output(model,model_data);
  sim_data.s = (v(:)') * sim_data.uh.dofs;
end;

model.decomp_mode = old_mode;