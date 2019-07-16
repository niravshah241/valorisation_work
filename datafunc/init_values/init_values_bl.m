function U0 = init_values_bl(glob,params)
%function U0 = init_values_bl(glob,params)
%
% function constructing the initial values of the convection diffusion
% problem in the specified global points glob and parameters.
% It returns an initial data function that is mosly homogeneous with several
% (seven) blob like structures of higher concentration that drops exponentially
% away from their centres.
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


% Martin Drohmann 23.9.2009

% glob column check
if params.debug
  if ~isempty(glob) && size(glob,1) < size(glob,2)
    warning('coordinates in variable glob are given row-wise, but expected them to be column-wise');
    if params.debug > 2
      keyboard;
    end
  end
end

decomp_mode = params.decomp_mode;

if ~isfield(params, 'c_width')
  params.c_width = 0.35;
end
if ~isfield(params, 'c_height')
  params.c_height = 0.50;
end

xoff = params.c_width/2;
yoff = params.c_height/2;

if decomp_mode == 2
  U0 = [params.c_high, params.c_low];
else
  X = glob(:,1);
  Y = glob(:,2);

  square = X>0.4-xoff & X<0.4+xoff & Y>0.5-yoff & Y < 0.5+yoff;

  if decomp_mode == 0
    U0 = zeros(length(X),1);
    U0(square) = params.c_high;
    U0(U0 == 0) = params.c_low;

  elseif decomp_mode == 1
    U0 = cell(2,1);
    U0{1} = zeros(length(X),1);
    U0{2} = zeros(length(X),1);
    U0{1}(square) = 1;
    U0{2} = 1 - U0{1};
  else
    error(['decomp_mode number ' decomp_mode, ' is unknown.']);
  end
end

