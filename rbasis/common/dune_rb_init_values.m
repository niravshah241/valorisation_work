function a0 = dune_rb_init_values( model, detailed_data )

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
  a0 = model.coeff_ops.u0_ptr(model.mu);
else
  if model.decomp_mode~=1
    model.mexptr('set_mu', model.mu);
  end
  a0 = model.mexptr('rb_init_values', model.decomp_mode);
end

end
%| \docupdate 
