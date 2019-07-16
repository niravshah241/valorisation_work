function Udirichlet = dirichlet_values_weighted_boxes(glob, params)
% function Udirichlet = dirichlet_values_weighted_boxes(glob, params)
%
% function computing dirichlet-values by pointwise evaluations  
%
% values constant params.c_dir in box given by
%       (beta=1):  params.dir_box_xrange{1} and dir_box_yrange{1} 
%       (beta=0):  params.dir_box_xrange{2} and dir_box_yrange{2} 
%            for intermediate betas, the two boxes are weighted
%
% required fields of params:
%   dir_middle : y coordinate where upper and lower regions are separated
%   c_dir      : dirichlet value in box given by beta and
%                'dir_box_{xrange,yrange}'
%   beta       : dirichlet values == params.c_dir in box given by
%       - '(beta=1)'  params.dir_box_xrange{1} and dir_box_yrange{1}
%       - '(beta=0)'  params.dir_box_xrange{2} and dir_box_yrange{2}
%       -     for intermediate betas, the two boxes are weighted
%   dir_box_xrange : cell array with two cells for length of the two boxes in
%                    x-direction
%   dir_box_yxange : cell array with two cells for length of the two boxes in
%                    y-direction

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


%
% in 'coefficients' mode, the parameters in brackets are empty

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
  if ~ismember('beta',params.mu_names) 
    Udirichlet = params.c_dir;
  else
    Udirichlet = params.c_dir * [params.beta; 1-params.beta];
  end;
else % decomp_mode < 2
  X = glob(:,1);
  Y = glob(:,2);

  if decomp_mode == 0
    Udirichlet = zeros(size(X)); 
    xrange = params.dir_box_xrange{1};
    yrange = params.dir_box_yrange{1};
    fi = X > xrange(1) & ...
      X < xrange(2) & ...
      Y > yrange(1) & ...
      Y < yrange(2);
    Udirichlet(fi) = params.c_dir*params.beta;
    xrange = params.dir_box_xrange{2};
    yrange = params.dir_box_yrange{2};
    fi = X > xrange(1) & ...
      X < xrange(2) & ...
      Y > yrange(1) & ...
      Y < yrange(2);
    Udirichlet(fi) = params.c_dir*(1-params.beta);
  elseif decomp_mode == 1
    % two components if beta is in mu, otherwise one component
    Udirichlet = cell(2,1);
    for q=1:2
      Udirichlet{q} = zeros(size(X)); 
      xrange = params.dir_box_xrange{q};
      yrange = params.dir_box_yrange{q};
      fi = X > xrange(1) & ...
        X < xrange(2) & ...
        Y > yrange(1) & ...
        Y < yrange(2);
      Udirichlet{q}(fi) = 1;
    end;
    if ~ismember('beta',params.mu_names) 
      % merge to one component
      Udirichlet = {params.beta*Udirichlet{1} + ...
        (1-params.beta)*Udirichlet{2}};
    end;
  else
    error('decomp_mode unknown');
  end
end
%| \docupdate 
