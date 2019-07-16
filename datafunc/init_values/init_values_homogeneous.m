function U0 = init_values_homogeneous(glob,params)
%function U0 = init_values_homogeneous(glob,params)
%
% function constructing the initial values of the convection diffusion
% problem in the specified global points glob and parameters.
% It returns a constant function scaled by a scalar params.c_init.
%
% required fields in params
%    c_init:   constant for homogeneous initial data to be returned
%
% in 'coefficient' mode the model_data structure is empty

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


% Bernard Haasdonk 22.7.2006

% glob column check
if params.debug
  if ~isempty(glob) && size(glob,1) < size(glob,2)
    warning('coordinates in variable glob are given row-wise, but expected them to be column-wise');
    if params.debug > 2
      keyboard;
    end
  end
end
if nargin ~= 2
  error('wrong number of parameters!');
end;

% affine decomposition is supported by init_data, so automatically
% satisfied.

decomp_mode = params.decomp_mode;

if decomp_mode     == 2
  U0 = params.c_init;
elseif decomp_mode   == 0
  U0 = params.c_init * ones(length(glob),1);
elseif decomp_mode == 1
  U0 = { ones(length(glob),1) };
else
  error(['decomp_mode number ', params.decomp_mode, ' is unknown.']);
end;
%| \docupdate 
