function diffusivity = diffusivity_richards_nonlinear(glob, params)
% function diffusivity = diffusivity_richards_nonlinear(glob, params)
% function computing the nonlinear diffusivity tensor of a richards problem by
% pointwise evaluation in the point sequences indicated by global coordinates
% in the columns of the matrix glob.
%
% It returns a gradient tensor for the Richards equation.
%
% required fields of params:
%      k                      : scalar scaling factor for diffusivity
%      gravity                : scalar modelling the physical gravity
%      richards_perm_ptr      : function pointer for permeability function
%      richards_retention_ptr : function pointer evaluating the retention curve
%      U                      : DoF vector of discrete solution
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
X=glob(:,1);
Y=glob(:,2);
U = params.U;
%vlen = size(U,1);
p_mu = spline_select(params);
U = U - params.gravity * Y.*(1+ppval(p_mu, X));
Kaval = params.k * real(params.richards_perm_ptr(U)) .* real(params.richards_retention_ptr(U));
diffusivity.K = Kaval;
diffusivity.epsilon = max(Kaval);
%    Utemp = U' * params.k + 0.0004;
%    diffusivity.K = spdiags(reshape([Utemp;zeros(1,vlen)],2*vlen,1),0,2*vlen,2*vlen);
%    diffusivity.epsilon = max(Utemp);

decomp_mode = params.decomp_mode;
if decomp_mode>0 && respected_decomp_mode==0
  error('function is nonlinear and does not support affine decomposition!');
end;

