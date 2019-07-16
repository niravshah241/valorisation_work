function U0 = init_values_wave(glob,params)
%function U0 = init_values_wave(glob,params)
%
% function constructing the initial values of the convection diffusion
% problem in the specified global points glob and parameters.
% It constructs a wave of values between 0 and c_init 
%           c_init* 0.5 * (sin(freq_x * x + freq_y * y) + 1)
%
% required fields in params:
% c_init            : maximum value in field
% freq_x            : frequency in x-direction
% freq_y            : frequency in y-direction
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
decomp_mode = params.decomp_mode;

if decomp_mode == 0
  U0 = params.c_init * 0.5 * (sin(params.freq_x * glob(:,1) + ...
			    params.freq_y * glob(:,2)  ) + 1 );
elseif decomp_mode == 1 
  % single component independent of mu_names
  if ismember('freq_x',params.mu_names) || ...
  ismember('freq.y',params.mu_names)
    error('affine decomp with respect to mu_names not possible!');
  end;
  U0 = {0.5 * (sin(params.freq_x * glob(:,1) + ...
	     params.freq_y * glob(:,2)  ) + 1 )};      
elseif  decomp_mode == 2 
  U0 = params.c_init;
else
  error(['decomp_mode number ', params.decomp_mode, ' is unknown.']);
end;

%| \docupdate 
