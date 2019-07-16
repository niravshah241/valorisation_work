function U = dirichlet_values_weighted_boxes_components(glob, params)
%function U = dirichlet_values_weighted_boxes_components(glob, params)
%
% function evaluating a function in the list of global coordinates
% specified in the columns of glob. Result is a cell array of
% matrices of dirichlet value components results as columns.
%
% values constant params.c_dir in box given by
%       (beta=1):  model.dir_box_xrange{1} and dir_box_yrange{1} 
%       (beta=0):  model.dir_box_xrange{2} and dir_box_yrange{2} 
%            for intermediate betas, the two boxes are weighted

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


% Bernard Haasdonk 27.8.2009

% glob column check
if params.debug
  if ~isempty(glob) && size(glob,1) < size(glob,2)
    warning('coordinates in variable glob are given row-wise, but expected them to be column-wise');
    if params.debug > 2
      keyboard;
    end
  end
end
U = cell(2,1);
for q=1:2
  U{q} = zeros(size(glob,1)); 
  xrange = params.dir_box_xrange{q};
  yrange = params.dir_box_yrange{q};
  fi = (X > xrange(1) & ...
        X < xrange(2) & ...
        Y > yrange(1) & ...
        Y < yrange(2));
  U{q}(fi) = 1;
end;
if ~ismember('beta',params.mu_names) 
  % merge to one component
  U = {params.beta*U{1} + ...
       (1-params.beta)*U{2}};
end;%| \docupdate 
