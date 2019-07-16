function detailed_data = basisgen_fixed(model,detailed_data)
%function detailed_data = basisgen_fixed(model,detailed_data)
% function generating a reduced basis with a fixed training set approach.
%
% The RB is started with initial data by varying all parameters over
% 'detailed_data.RB_info.M_train'.  The basis is extended by searching for the
% vector in the training set with the maximum posterior error estimator or true
% error. For this parameter a single vector is chosen for basisextension.  The
% resulting reduced basis vectors are returned in RB which is a matrix of
% columnwise reduced basis vector-DOFS. The algorithm supports restart, i.e.
% continuing of prescribed initial basis RB.
%
% generated fields of detailed_data:
%           RB : collection of orthonormal reduced basis DOF
%                vectors
%           RB_info: datastructure with basis-generation
%                information, see RB_greedy_extension() for details
%
% required fields of detailed_data:
%           RB_info.M_train : matrix with columnwise parameter
%                             vectors to use in greedy search
%           grid : numerical grid, is constructed if not existent.
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


% Bernard Haasdonk 27.3.2007

% We assume, that detailed_data is reasonably filled, e.g. contains
% the grid.
%
% generate grid
%if ~isfield(detailed_data,'grid') || isempty(detailed_data.grid)
%  detailed_data.grid = construct_grid(model); 
%end;

%%%%%%%%%%%%%%% specific settings for fixed algorithm

if (model.get_rb_size(model,detailed_data)==0)  
  % generate RBstart: all initial data constellations
  RB = model.rb_init_data_basis(model, detailed_data);
  detailed_data = model.set_rb_in_detailed_data(detailed_data,RB);
%  if (get_rb_size(detailed_data)==0) % check if reduced basis is empty
%  detailed_data = model.RB_init_data_basis(model, detailed_data, ...
%                                           detailed_data.RB_info.M_train);
end;
%prepare model 
%model.Nmax = size(detailed_data.RB,2);

disp('Starting RB extension loop');

% start time measurement
tic;
start_time_basisgen_fixed = tic;

% extend basis 
detailed_data = RB_greedy_extension(model, ...
                                    detailed_data);

disp(['Generated RB basis on fixed M_train with ',...
      num2str(model.get_rb_size(model,detailed_data)), ...
      ' basis vectors.']);
t = toc(start_time_basisgen_fixed);
disp(['Runtime = ',num2str(t)]);



