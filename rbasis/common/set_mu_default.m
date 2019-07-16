function model = set_mu_default(model,mu,dummy)
%function model = set_mu_default(model,mu)
% 
% function setting the parameter vector mu in the model struct
% required fields of model:
% mu_names : cell array of strings, indicating the fields, which are set
%            by the current routine
%
% parameters:
% mu : new parameter vector `\mu`

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

if nargin == 2
  dummy = 0;
end

if nargout == 0
  warning('set_mu called without return argument!');
end

if length(mu(:))~=length(model.mu_names)
  error('dimensionality of mu does not fit to model.mu_names');
end;

for i=1:length(mu(:))
  model.(model.mu_names{i}) = mu(i);
end;

