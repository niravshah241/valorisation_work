function f = output_function_box_mean(glob,params)
%function f = output_function_box_mean(glob,params)
%
% function computing the integral kernel of a linear output
% functional given as integral:
%
%    output(u) = int_Omega f u dx
% 
% the function f is represented by this function, the integration
% performed in suitable discretization routines later, e.g. 
% fv_operators_output. glob is assumed to be a matrix of columnwise
% global coordinate vectors. Res is a matrix with columnwise
% results of application of f on all columns of glob.
%
% box_mean:  for computation of mean(U) over box specified by
% params, i.e. f indicator function of the box scaled with the
% inverse volume of the box.
%
% note, that this functional is only exact, if the boundary of the
% box is exaclty corresponding ot element boundarys of the
% numerical grid.
%
% required fields of params:
%   sbox_xmin, sbox_xmax sbox_ymin, sbox_ymax: coordinates of the
%   box over which averaging is performed. 
% 
% These parameters should not be under parameter variation control,
% but be constant throughout the parameter variation.
% glob column check

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

if params.debug
  if ~isempty(glob) && size(glob,1) < size(glob,2)
    warning('coordinates in variable glob are given row-wise, but expected them to be column-wise');
    if params.debug > 2
      keyboard;
    end
  end
end

% Bernard Haasdonk 16.5.2008

decomp_mode = params.decomp_mode;

if decomp_mode == 2
  f = 1; % one component weight 1, i.e. no parametrization
else
  
  boxvolinv = 1/((params.sbox_xmax-params.sbox_xmin)*...
		 (params.sbox_ymax-params.sbox_ymin));
  
  I =  (glob(:,1)<=params.sbox_xmax) & (glob(:,1)>=params.sbox_xmin) ...
       & (glob(:,2)>=params.sbox_ymin) & (glob(:,2)<=params.sbox_ymax);

%  I2 = (glob(:,1)>=params.sbox_xmin) & (glob(:,1)<=params.sbox_xmax) & ...
%       (glob(:,2)<=params.sbox_ymax) & (glob(:,2)>=params.sbox_ymin);

%  I3 = (glob(:,1)<=params.sbox_xmax) & (glob(:,1)>=params.sbox_xmin) ...
%       & (glob(:,2)<=params.sbox_ymax) & (glob(:,2)>=params.sbox_ymin);

%  keyboard;
  
  % function simply is 0 if glob is not in the box, 
  % 1/vol(box) if glob is in the box.
  f = I.*boxvolinv;  
  
  if decomp_mode == 1  % single component in cell array
    f = {f};
  else % decomp_mode == 0
    %f = f;
  end;
end;



%| \docupdate 


