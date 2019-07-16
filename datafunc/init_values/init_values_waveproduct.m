function U0 = init_values_waveproduct(glob, params)
% function U0 = init_values_waveproduct(glob, params)
% product of two sinus waves in both coordinate directions
%
% computes a product of axis-dependent waves of values between 'c_init_min'
% and 'c_init_max'
%
% required fields of params:
% c_init_max     : maximum value in field
% c_init_min     : minimum value in field
% c_init_freq_x  : frequency in x-direction
% c_init_freq_y  : frequency in y-direction
% c_init_phase_x : phase shift
% c_init_phase_y : phase shift

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

%  if ismember('c_init_freq_x',params.mu_names) | ...
%       ismember('c_init_freq_y',params.mu_names) | ...
%       ismember('c_init_min',params.mu_names) | ...
%       ismember('c_init_max',params.mu_names) | ...
%       ismember('c_init_phase_x',params.mu_names) | ...
%       ismember('c_init_phase_y',params.mu_names) 
%    error('affine decomp with respect to mu_names not possible!');
%  end;

if decomp_mode == 2
  U0 = 1;
else
  X = glob(:,1);
  Y = glob(:,2);

  Uinit = params.c_init_min + ...
          (params.c_init_max-params.c_init_min) * ...
          0.5 * (sin_sym(params.c_init_freq_x * X(:) +  params.c_init_phase_x) ...
                 .* sin_sym(params.c_init_freq_y * Y(:) +  params.c_init_phase_y) ...
                 + 1);
  if decomp_mode == 0
    U0 = Uinit;
  elseif decomp_mode == 1
    U0 = {Uinit};
  end
end
%| \docupdate 
