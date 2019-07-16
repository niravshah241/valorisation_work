function U0 = init_values_transformed_blobs_richards(glob,params)
%function U0 = init_values_transformed_blobs_richards(glob,params)
%
% function constructing the initial values of the convection diffusion
% problem in the specified global points glob and parameters.
% It returns an initial data function that is mosly homogeneous with several
% (seven) blob like structures of higher concentration that drops exponentially
% away from their centres.
%
% required fields in params
%    c_init:   constant for homogeneous initial data to be returned
%    blob_height: addend to c_init for the maximum concentration of the blobs
%    gravity:     physical gravity
%    hill_height: geometry parameter for hill of upper gemeotries' boundary
%
% in 'coefficient' mode the model_data structure is empty

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
decomp_mode = params.decomp_mode;

if any(ismember(fieldnames(params), 'blob_height'))
  h = params.blob_height;
else
  h = 0.15;
end

if decomp_mode == 2
  U0 = [ h, h, h, h, h, h, h, params.c_init];
  U0 = [ U0, params.gravity, params.gravity * params.hill_height ];
else
  X = glob(:,1);
  Y = glob(:,2);
  B{1} = ((X(:) - 0.5).^2 +(Y(:)-0.7).^2);
  B{2} = ((X(:) - 0.8).^2 +(Y(:)-0.9).^2);
  B{3} = ((X(:) - 0.8).^2 +(Y(:)-0.2).^2);
  B{4} = ((X(:) - 0.2).^2 +(Y(:)-0.2).^2);
  B{5} = ((X(:) - 0.2).^2 +(Y(:)-0.9).^2);
  B{6} = ((X(:) - 0.4).^2 +(Y(:)-0.4).^2);
  B{7} = ((X(:) - 0.6).^2 +(Y(:)-0.6).^2);
                  
  I=cell(1,7);
  % p(mu) = mu*(1-x)
  % p_mu = spline_select(params);
  for i=1:7
      I{i} = (B{i} < 0.05.^2);
  end
  if decomp_mode == 0
    U0 = h .* I{1} + params.c_init * ones(length(X),1);
    % TODO: replace ztrans with something derivated from myspline

    ztrans = Y.*(1-X);
    U0 = U0 + params.gravity * (Y + params.hill_height*ztrans);

    for i = 2:7
      U0 = U0 + h .* exp(-B{i}) .* I{i};
    end
  elseif decomp_mode == 1
    ztrans = Y.*(1-X);
    U0 = cell(1,8);
    U0{1} = I{1};
      for i = 2:7
      U0{i} = I{i};
      end
    U0{8} = ones(length(X),1);
    U0{9} = Y;
    U0{10} = ztrans;
  else
    error(['decomp_mode number ' params.decomp_mode, ' is unknown.']);
  end
end
%| \docupdate 
