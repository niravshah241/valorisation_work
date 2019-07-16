function RBinit = RB_init_data_basis(model,detailed_data)
%function RBinit = RB_init_data_basis(model,detailed_data)
%
% function generating an initial reduced basis by varying all mu-values
% in the columns of M and collecting the init-states. 
% detailed_data is assumed to contain the list of parameters in
% detailed_data.RB_info.M_train
%
% required fields of model:
%     init_values_algorithm : name of function computing the
%                     initial values
%     inner_product_matrix_algorithm : function giving the inner-product matrix

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

%grid = detailed_data.grid;

% generate RBinit: all initial data constellations
RBinit = [];
M = detailed_data.RB_info.M_train;
model.decomp_mode = 0;

for i = 1:size(M,2)
  model = model.set_mu(model,M(:,i));
  % use detailed_data as model_data here:
  U0 = model.init_values_algorithm(model,detailed_data);
  RBinit = [RBinit, U0];
end;

RBinit = model.orthonormalize(model, detailed_data, RBinit);
%  RBinit = delzerocolumns(RBinit);

% start with very simple RBinit if none was found yet
if(size(RBinit, 2) == 0)
  RBinit = ones(size(RBinit,1), 1);
end
disp(['found ',num2str(size(RBinit,2)),' basis functions for',...
      ' init data variation.']);
 
