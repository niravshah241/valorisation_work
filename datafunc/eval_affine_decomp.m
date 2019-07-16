function R = eval_affine_decomp(R_comp_ptr,R_coeff_ptr,varargin)
%function R = eval_affine_decomp(R_comp_ptr,R_coeff_ptr,varargin)
%
% function returning either the components, coefficients or their
% linear combination depending on params.decomp_mode, where
% params.decomp_mode is assumed to be the last argument in
% varargin. Both function pointers are assumed to have the
% identical command line syntax R_comp_ptr(varargin(:));

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


params = varargin{end};
if params.decomp_mode == 2
  R = R_coeff_ptr(varargin{:});
else
  if params.decomp_mode == 1
    R = R_comp_ptr(varargin{:});    
  else
    Rcomp = R_comp_ptr(varargin{:});    
    Rcoeff = R_coeff_ptr(varargin{:});
    R = lincomb_sequence(Rcomp,Rcoeff);
  end;
end;

%| \docupdate 
