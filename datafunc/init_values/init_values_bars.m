function U0 = init_values_bars(glob,params)
%function U0 = init_values_bars(glob,params)
%
% function constructing the initial values of the convection diffusion
% problem in the specified global points glob and parameters.
% It returns an initial data function that is mosly homogeneous with several
% (seven) blob like structures of higher concentration that drops exponentially
% away from their centres.
%
% required fields in params
%    c_init:   constant for homogeneous initial data to be returned
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

if ~isfield(params, 'bars_width')
  params.bars_width = 0.10;
end

if decomp_mode == 2
  U0 = [params.c_init, [1, 1/2, 1/3, 1/4, 1/5]*0.5];
else
  X = glob(:,1);
  Y = glob(:,2);
  offset1 = 1.1-params.bars_width/2;
  offset2 = 1.1+params.bars_width/2;
  if decomp_mode == 0
    U0 = params.c_init * ones(length(X),1);
    Uleft = zeros(size(X));
    for i = 1:3
%      U0 = U0 + 0.5/i*(X(:) > -0.123+0.2*i & X(:) < -0.077+0.2*i & Y(:) < 0.5);
      Uleft = Uleft + 0.5/i*(X > offset1-0.2*i & X < offset2-0.2*i & Y > 0.5 & X > 0.5) ...
                    + 0.5/(6-i)*(X > offset1-0.2*i & X < offset2-0.2*i & Y < 0.5 & X > 0.5);
%      U0 = U0 + 0.5/i*(X(:) > 1.075-0.2*i & X(:) < 1.125-0.2*i & Y(:) > 0.5);
%      U0 = U0 + 0.5/(6-i)*(X(:) > 1.075-0.2*i & X(:) < 1.125-0.2*i & Y(:) < 0.5);
    end
    ysize = size(Uleft, 2);
    xsize = size(Uleft, 1);

    Uright(xsize:-1:1, ysize:-1:1) = Uleft;
    U0 = U0 + Uleft(:) + Uright(:) + 0.5/3 * (X==0.5);
  elseif decomp_mode == 1
    U0 = cell(7,1);
    U0{1} = ones(length(X),1);
    Uleft = cell(6,1);
    Uright = cell(6,1);
    for i = 1:6
      Uleft{i} = zeros(size(X));
    end
    for i = 1:3
      Uleft{i} = Uleft{i} + ...
        (X > offset1-0.2*i & X < offset2-0.2*i & Y > 0.5 & X > 0.5);
      Uleft{6-i} = Uleft{6-i} + ...
        (X > offset1-0.2*i & X < offset2-0.2*i & Y < 0.5 & X > 0.5);
    end
    ysize = size(Uleft{1}, 2);
    xsize = size(Uleft{1}, 1);
    for i = 1:6
      Uright{i}(xsize:-1:1, ysize:-1:1) = Uleft{i};
    end
    for i = 1:6
      U0{i+1} = Uleft{i}(:) + Uright{i}(:);
    end
    U0{4} = U0{4} + 0.5/3 * (X==0.5);
  else
    error(['decomp_mode number ' decomp_mode, ' is unknown.']);
  end
end

%| \docupdate 
