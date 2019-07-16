function [vel,lambda] = velocity_parabola(glob, params)
% function [vel,lambda] = velocity_parabola(glob, params)
% parabolic velocity field
%
% required fields of params:
%  c : velocity of parabola-profile in 'params.yrange',
%      constant in 'x'-direction

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

% determine affine_decomposition_mode as integer
decomp_mode = params.decomp_mode;

if decomp_mode==2
  vel = params.c;
  %| \todo lambda needs to be set to something more
  %reasonable for reduced simulations.
  lambda = 0;
else
  X = glob(:,1);
  Y = glob(:,2);
  if decomp_mode == 1
    vel = { [ zeros(length(X),1), Y(:) .* (1-Y(:)) ] };
    %| \todo lambda needs to be set to something more
    %reasonable for reduced simulations.
    lambda = 0;
  elseif decomp_mode == 0 % decomp_mode 0
    % determine lambda_jl to be globally constant such that
    %   lambda_jl * sup_u n_jl * f'(u) <= 1
    %  f'(u) = v(x,y)
    %  e.g. lambda := 1/sup|v(x,y)|
    lambda = 1/(abs(params.c)/4.0+ 1e-10); % some eps for divby0
    vel = [ zeros(length(X),1), params.c * Y(:) .* (1-Y(:)) ];
  else
    error('unknown decmposition mode!!');
  end
end

%| \docupdate 
