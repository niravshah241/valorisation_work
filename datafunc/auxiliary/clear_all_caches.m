function clear_all_caches
% function clear_all_caches
% This function clears the caches generated for caching of gradient data in
% evolution schemes.
%
% In detail the caches in inv_geo_trans_derivative() and diffusivity_cached()
% are cleared.

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

  inv_geo_trans_derivative();
  diffusivity_cached();
end

