function rb_sim_data = lin_stat_rb_reconstruction(model,detailed_data,rb_sim_data)
%function rb_sim_data = lin_stat_rb_reconstruction(model,detailed_data,rb_sim_data)
%
% (trivial) function computing a detailed reconstruction by linear
% combination of the coefficients in the simulation data with the
% orthonormal reduced basis RB

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


% Bernard Haasdonk 22.2.2011

if ~isfield(rb_sim_data,'uh')
  rb_sim_data.uh = femdiscfunc([],detailed_data.df_info);
end;  
rb_sim_data.uh.dofs = ...
    detailed_data.RB(:,1:length(rb_sim_data.uN)) * rb_sim_data.uN;
%| \docupdate 
