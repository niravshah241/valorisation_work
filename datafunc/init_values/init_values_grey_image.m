function U0 = init_values_grey_image(glob,params)
%function U0 = init_values_grey_image(glob,params)
%
% function constructing the initial values of the convection diffusion
% problem in the specified global points glob and parameters.
% load grey value image as initial data,
% a grey value image is loaded mapping 0 to 0 and 255 to c_init
% with bilinear interpolation 
%
% required fields in params:
% c_init             : parameter weighting the components
% c_init_filename    : name of image file
% c_init_xmin        : image is mapped to specified rectangle 
% c_init_xmax        : image is mapped to specified rectangle 
% c_init_ymin        : image is mapped to specified rectangle 
% c_init_ymax        : image is mapped to specified rectangle 
% c_init_xdivisions  : vector of divisions, i.e intervals,
%                      e.g. borders of letters in an image
%                      determines number of components. c_init = 0
%                      => full weight to leftmost part of image
%                      c_init = 1 => full weight to complete image
%
% in 'coefficient' mode the params structure is empty

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






ncomp = 1+length(params.c_init_xdivisions);
if decomp_mode ~= 2
  % determine components and perform linear combination
  img = imread(params.c_init_filename)/255;
  %erstmal ohne bilineare Interpolation:
  XI = (glob(:,1)-params.c_init_xmin)/(params.c_init_xmax- ...
				  params.c_init_xmin)* size(img,2);
  XI = max(XI,0);
  XI = min(XI,size(img,2));
  XI = ceil(XI);
  fi = find(XI==0);
  if ~isempty(fi)
    XI(fi)= 1;
  end;    
  YI = (glob(:,2)-params.c_init_ymin)/(params.c_init_ymax- ...
				  params.c_init_ymin)* size(img,1);
  YI = max(YI,0);
  YI = min(YI,size(img,1));
  % reflect y coordinate
  YI = size(img,1)-YI;
  YI = ceil(YI);
  fi = find(YI==0);
  if ~isempty(fi)
    YI(fi)= 1;
  end;
  
  ind = sub2ind(size(img),YI,XI);
  I = img(ind);
  
  U0comp = cell(ncomp,1);
  xdivs = [0;params.c_init_xdivisions(:);params.c_init_xmax];
  for i = 1:ncomp
    U0 = zeros(length(glob(:,1)),1);
    fi = find( (glob(:,1)<xdivs(i+1)) & (glob(:,1)>=xdivs(i)));
    if ~isempty(fi)
      U0(fi) = I(fi);
    end;
    U0comp{i} = U0;
  end;
end;

% determine coefficients:
coeff = zeros(ncomp,1);
coeff = (0:(ncomp-1))/ncomp;
coeff = (params.c_init - coeff)*ncomp;
coeff = min(coeff,1);
coeff = max(coeff,0);

% now return the correct quantity
if decomp_mode == 2 
  U0 = coeff;
elseif decomp_mode == 0 
  U0 = lincomb_sequence(U0comp,coeff);
elseif decomp_mode == 1 
  U0 = U0comp;
else
  error(['decomp_mode number ', params.decomp_mode, ' is unknown.']);
end;
respected_decomp_mode = 1;

