function [sim_data,tictoc] = load_detailed_simulation(m,savepath,params)
% function [sim_data,tictoc] = load_detailed_simulations(m,savepath,params)
% load single trajectory of previously saved results.
%
% If a routine uses this method, a previous save_detailed_simulations() should
% have been performed to guarantee consistency of the saved data with the
% current grid and parameters.
%
% \note No consistency check is performed.
%
% Parameters:
%  m:        is the index of the function to be looked for.
%  savepath: is a path relative to 'RBMATLABTEMP', where the functions are
%            stored by save_detailed_simulations()
%  params:   is a control structure (currently unused)
%
% Return values:
%  sim_data:  is the loaded simulation data structure.
%  tictoc:    is the time in seconds it took to compute the detailed
%             simulation.

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


% Bernard Haasdonk 14.5.2007

sp = fullfile(rbmatlabtemp,savepath);
load(fullfile(sp,['detail',num2str(m),'.mat']));

if isempty(tictoc)
  tictoc = 0;
end

