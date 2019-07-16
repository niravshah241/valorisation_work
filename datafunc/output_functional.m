function s = output_functional(model,model_data,U)
%function s = output_functional(model,model_data,U)
%
% function computing an output functional from the discrete
% function U. Grid is needed to have the space discretization 
% information about U.
% return values are s.value and s.l2norm the induced functional
% l2-norm of the output functional.

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


% 'box_mean':  s = mean(U) over box specified by 
%
% required fields of model:
% name_output_functional  : 'box_mean'
% name_discfunc_element_mean : name of function, that computes the
%                  mean values of a discrete function on a given
%                  set of grid elements 
%
% if name_init_values == 'box_mean'
%   sbox_xmin, sbox_xmax sbox_ymin, sbox_ymax: coordinates of the
%   box over which averaging is performed. 
% 
% The implementation is a bit inprecise: The function U is only
% evaluated in the cell-centroids and weighted by the cell
% volume. This is summed for all cells with centroid in the specified box... 
%
% these parameters should not be under parameter variation control,
% but be constant throughout the parameter variation.

% Bernard Haasdonk 16.5.2008

disp('deprecated, please use pointer to suitable output_functional_*')

grid = model_data.grid;

if isequal(model.name_output_functional,'box_mean')  
  % definition of the box
  xmin = model.sbox_xmin;
  xmax = model.sbox_xmax;
  ymin = model.sbox_ymin;
  ymax = model.sbox_ymax;
  
  I = find((grid.CX<=xmax) & (grid.CX>=ymin) & (grid.CY>=ymin) ...
	   & (grid.CY<=ymax));
  
  U_means = model.name_discfunc_element_mean(model,model_data,U,I); 
  
  Apart = sum(grid.A(I));
  s.value = sum(U_means.*grid.A(I))/Apart;
  
  % l2norm of output functional: sup_|u|=1 s(u), i.e.
  % the discrete function u that has constant value c in the sensor
  % domain => l2norm(u)=sqrt(Apart * c^2) == 1
  % hence c = sqrt(1/Apart) and  s(u) = ||s|| = c  
  s.l2norm = sqrt(1/Apart);
    
else
  error('unknown name_output_functional');
end;

%| \docupdate 
