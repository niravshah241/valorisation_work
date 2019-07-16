function l2error = dune_lin_evol_rb_recons_and_compare(model,rb_sim_data)
%function dune_lin_evol_rb_recons_and_compare(model, rb_sim_data)
%
% reconstruct the reduced function and compare it to a full solution
%
% Sven Kaulmann 18.12.2009
%

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

  
  % second argument is detailed data but isn't needed
 l2error = model.mexptr('rb_reconstruct_and_compare', rb_sim_data, rb_sim_data);
end