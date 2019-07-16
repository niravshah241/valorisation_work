function diffusivity = diffusivity_homogeneous(glob, params)
% function diffusivity = diffusivity_homogeneous(glob, params)
%
% function computing the diffusivity pointwise evaluation in the point
% sequences indicated by global coordinates in the columns of the matrix glob.
% It returns a homogeneous diffusion coefficient.
%
% required fields of params:
%       k : scalar scaling factor
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
decomp_mode = params.decomp_mode;
diffusivity = [];
diffusivity.epsilon = 0;

if decomp_mode == 2 % decomp_mode ==2, single component
  diffusivity = params.k;
elseif decomp_mode == 0
  diffusivity.epsilon = params.k;
  diffusivity.K = params.k * ones(length(glob),1);
elseif decomp_mode == 1
  d = [];
  % single component independent whether k in mu
  d.epsilon = params.k;
  d.K = ones(length(glob),1);
  diffusivity = {d};
else
  error('unknown decomp mode');
end;
%| \docupdate 
