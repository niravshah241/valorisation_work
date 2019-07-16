function model = dune_set_mu(model, mu, skip_send_to_dune)
%function dune_set_mu(model, mu)
% sets the parameter mu in dune-rb

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


model.mu = mu;

if nargout == 0
  warning('dune_set_mu called without return argument!');
end

if nargin==2 || ~skip_send_to_dune
  model.mexptr('set_mu', mu);
end


end%| \docupdate 
