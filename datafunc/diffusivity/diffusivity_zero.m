function diffusivity = diffusivity_zero(glob, params)
% function diffusivity = diffusivity_zero(glob, params)
%
% function computing the diffusivity pointwise evaluation in the point
% sequences indicated by global coordinates in the columns of the matrix glob.
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
diffusivity.K       = 0;
diffusivity.epsilon = 0;

%| \docupdate 
