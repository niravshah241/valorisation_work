function U0 = init_values_blobs(glob,params)
%function U0 = init_values_blobs(glob,params)
%
% function constructing the initial values of the convection diffusion
% problem in the specified global points glob and parameters.
% u_0 = beta*exp(-gamma*((x-0.25).^2+(y-0.5).^2))* Xi(B_r) ...
%                  + (1-beta)*exp(-gamma*((x-0.25).^2+(y-0.7).^2))* Xi(B_r)
% required fields:
% radius             : radius of blob-cutoff circle
% gamma              : spread of initial-value gaussian
% beta               : value between 0 and 1 weighting two gauss-blobs
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

if decomp_mode == 2
  U0 = [params.beta, 1-params.beta];
else
  X = glob(:,1);
  Y = glob(:,2);

  B1 = ((X(:)-0.25).^2+(Y(:)-0.5).^2);
  B2 = ((X(:)-0.25).^2+(Y(:)-0.7).^2);
  I1 = (B1<params.radius.^2);
  I2 = (B2<params.radius.^2);
  if decomp_mode == 0 
    U0 = params.beta*exp(-params.gamma*B1).*I1 ...
      + (1-params.beta)* exp(-params.gamma*B2).*I2;
  elseif decomp_mode == 1 
    if ismember('gamma',params.mu_names) || ...
        ismember('radius',params.mu_names)
      error('affine decomp with respect to mu_names not possible!');
    end;
    % two components
    U0{1} = exp(-params.gamma*B1).*I1;
    U0{2} = exp(-params.gamma*B2).*I2;
  end
end
%| \docupdate 
