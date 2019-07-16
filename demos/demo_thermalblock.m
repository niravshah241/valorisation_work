function demo_thermalblock
% function demo_thermalblock
%
% rb demo gui with the termal block. currently simple
% snapshot-basis is used, so no global good approximation, but
% functionality and error estimation can nicely be demonstrated.

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


% B.Haasdonk, I. Maier 26.04.2011

disp('demo_rb_gui:')  
model = thermalblock_model;
model_data = gen_model_data(model);
detailed_data = gen_detailed_data(model,model_data);
plot_params.axis_tight = 1;
plot_params.yscale_uicontrols = 0.5;
demo_rb_gui(model,detailed_data,[],plot_params);
  