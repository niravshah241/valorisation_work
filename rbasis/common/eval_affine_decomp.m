function R = eval_affine_decomp(R_comp_ptr,R_coeff_ptr,model,model_data)
%function R = eval_affine_decomp(R_comp_ptr,R_coeff_ptr,model,model_data)
%
% function returning either the components, coefficients or their
% linear combination depending on model.decomp_mode, where
% expected argument list: R_comp_ptr(model,model_data);
% R_coeff_ptr(model);

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


if model.decomp_mode == 2
  R = R_coeff_ptr(model);
else
  if model.decomp_mode == 1
    R = R_comp_ptr(model,model_data);    
  else
    Rcomp = R_comp_ptr(model,model_data);    
    Rcoeff = R_coeff_ptr(model);
    R = lincomb_sequence(Rcomp,Rcoeff);
  end;
end;

%| \docupdate 
