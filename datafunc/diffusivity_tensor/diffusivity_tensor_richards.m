function diff = diffusivity_tensor_richards(glob,U,params,callerid)
%function diffusivity = diffusivity_tensor_richards(glob,U,params,callerid)
% function computing the diffucivity tensor of a diffusion problem by pointwise
% evaluation in the point sequences indicated by 'glob' and dof vector 'U'.
%
% Parameters:
%   U: Dof vector which is unneeded here, as the tensor is linear.
%   params: structure controlling the parameters of this function
%   callerid: The callerid needs to be set in case the diffusivity tensor() is
%             cached for later time steps. It adds a unique item to the
%             argument list of a later inv_geo_trans_derivative() call. That
%             allows correct hashing.
%
% generated fields of diff:
%     epsilon: upper bound on diffusivity value
%     K      : diffusivity tensor, i.e. a sparse square matrix of size 
%              2*dim x 2*dim
%

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


% Bernard Haasdonk 11.4.2006

% determine affine_decomposition_mode as integer  
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
% flag indicating whether the computation respected the decomposition

if decomp_mode~=0 && ~isempty(intersect(params.mu_names,'hill_height'))
  error('affine decomposition not implemented!!');
end;

diff = [];
diff.epsilon = 0;

%  [res1, res2] = inv_geo_trans_derivative(glob,{(1), (2)},{(1), (2)},params);
%  row1 = [res1{1}, res1{2}];
%  row2 = [res2{1}, res2{2}];

%  vlen = size(row1,1);
%  temp0  = reshape([ sum(row1 .* row1, 2), sum(row2 .* row2, 2) ]',2*vlen,1);
%  tempm1 = reshape([ sum(row2 .* row1, 2), zeros(vlen,1) ]', 2*vlen, 1);
%  tempp1 = reshape([ zeros(vlen,1), sum(row1 .* row2, 2) ]', 2*vlen, 1);
%  temp1 = [ sum(row1 .* row1, 2), sum(row1 .* row2, 2) ];
%  temp2 = [ sum(row2 .* row1, 2), sum(row2 .* row2, 2) ];
%  tempD  = spdiags([tempm1,temp0,tempp1],-1:1,2*vlen,2*vlen);
tempD = real(diffusivity_cached(glob,params,callerid));

%  min(min(abs(atemp1 - temp1) < 1e-16))
%  min(min(abs(atemp2 - temp2) < 1e-14))
%  keyboard;

%  d.K1   = repmat(diff_k.K,1,2) .* temp1;
%  d.K2   = repmat(diff_k.K,1,2) .* temp2;


%  d.K     = diff_k.K * tempD;
d.K     = tempD;
  
%  d.K1    = params.k * temp1;
%  d.K2    = params.k * temp2;
%  keyboard
d.epsilon = max(max(d.K));

diff = d;


