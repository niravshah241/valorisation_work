function U0 = init_values_gradient_box(glob,params)
%function U0 = init_values_gradient_box(glob,params)
% init value function with a gradient box profile in Y-direction as output
%
% function constructing the initial values of the convection diffusion
% problem in the specified global points glob and parameters.
% ``u_0 = ( c_0 (y-y_{\min}) + c_1 (y_{\max}-y) ) \cdot
% \frac{1}{y_{\max}-y_{\min}} \cdot \chi_{B}``
% where `B=[x_{\min},x_{\max}] \times [y_{\min},y_{\max}]` denotes the
% gradient box with `(x_{\min},x_{\max}) = (0.35,0.55)` and
% `(y_{\min},y_{\max}) = (0.2,0.9)`.
%
% required fields of params:
%  c_init_up  : constant value at upper boundary of gradient box `c_0`
%  c_init_low : constant value at lower boundary of gradient box `c_1`

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

decomp_mode = params.decomp_mode;

if decomp_mode == 2 
  U0 = [params.c_init_up, params.c_init_lo];
else
  width = 0.2;
  xmid = 0.45;
  ymid = 0.55;
  height = 0.7;
  ymin = ymid-height/2;
  ymax = ymid+ height/2;
  if decomp_mode == 0
    BX = abs(glob(:,1)-xmid);
    BY = abs(glob(:,2)-ymid);
    I = (BX<width/2) & (BY<height/2);
    % y = ymin =>  0 + c_init_lo
    % y = ymax =>  c_init_up + 0
    U0 = (params.c_init_up * (glob(:,2)-ymin) + ...
          + params.c_init_lo * (ymax-glob(:,2))) .* I / (ymax-ymin);
  elseif decomp_mode == 1
    % two components
    BX = abs(glob(:,1)-xmid);
    BY = abs(glob(:,2)-ymid);
    I = (BX<width/2) & (BY<height/2);
    U0{1} = (glob(:,2)-ymin) .* I / (ymax-ymin);
    U0{2} = (ymax-glob(:,2)) .* I / (ymax-ymin);
  end
end

