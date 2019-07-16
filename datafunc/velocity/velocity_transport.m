function [vel,lambda] = velocity_transport(glob, params)
%function [vel,lambda] = velocity_transport(glob, params)
%
% function evaluating a function in the list of global coordinates
% specified in the columns of glob. Result is a matrix of velocity
% vectors as columns of vel.
%
% Linear combination of components by coefficients then yields the
% complete evaluation.

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

lambda = 0;
decomp_mode = params.decomp_mode;

if decomp_mode == 2
  vel = 1; % single component factor 1
elseif decomp_mode == 1 % components:
  vel = cell(1,1);
  vel{1} = zeros(length(glob),2);
  vel{1}(:,1) = params.transport_x * ones(1,length(glob));
  vel{1}(:,2) = params.transport_y * ones(1,length(glob));
else % decomp_mode = 0, complete:
  vel = zeros(length(glob),2);
  vel(:,1) = params.transport_x * ones(1,length(glob));
  vel(:,2) = params.transport_y * ones(1,length(glob));
  lambda = max(max(vel))^-1; % correct up to a factor
end;

%| \docupdate
