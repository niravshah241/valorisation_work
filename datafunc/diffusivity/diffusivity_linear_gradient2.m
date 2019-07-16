function diffusivity = diffusivity_linear_gradient2(glob, params)
% function diffusivity = diffusivity_linear_gradient2(glob, params)
%
% function computing the diffusivity pointwise evaluation in the point
% sequences indicated by global coordinates in the columns of the matrix glob.
% It returns a gradient that increases linearly from model.diff_left to
% model.diff_right in X direction of the geometry.
%
% generated fields of diffusivity:
%     epsilon: upper bound on diffusivity value
%     K: vector with diffusivity values

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
decomp_mode = model.decomp_mode;
diffusivity = [];
diffusivity.epsilon = 0;

temp_m = (model.diff_right-model.diff_left) / ...
         (model.xrange(2) - model.xrange(1));
lin_grad = @(x) (x - model.xrange(1) );
if decomp_mode == 2
  diffusivity = temp_m;
elseif decomp_mode == 0
  diffusivity.epsilon = max(model.diff_left, model.diff_right);
  diffusivity.K = model.diff_left + (lin_grad(glob(:,1)) * temp_m);
elseif (decomp_mode == 1)
  d = [];
  d.epsilon = max(model.diff_left, model.diff_right);
  d.K = lin_grad(glob(:,1));
  diffusivity = {d};
else
  error('unknown decomp mode');
end
%| \docupdate 
