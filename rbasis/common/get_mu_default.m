function mu = get_mu_default(model)
%function mu = get_mu_default(model)
% 
% function getting the current parameter vector mu in the model struct
%
% required fields of model:
% mu_names : cell array of strings, indicating the fields, whcih are read
%            by the current routine

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

  
% Bernard Haasdonk 20.7.2006

mu = zeros(length(model.mu_names),1);
for i=1:length(mu(:))
  mu(i) = model.(model.mu_names{i});
end;

