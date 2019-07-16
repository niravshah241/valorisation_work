function model=model_default(model, T, nt)
% model = model_default(model)

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


model.t               = 0;
model.tstep           = 1;
model.decomp_mode     = 0;
model.verbose         = 0;
model.debug           = 0;
model.orthonormalize  = @model_orthonormalize_qr;
model.dt              = model.T / model.nt;
model.ei_time_indices = 1:model.nt+1;
model.mu              = zeros(size(model.mu_names));
%| \docupdate 
