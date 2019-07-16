function [enlarged, dummy] = dune_RB_extension_PCA_fixspace( model, detailed_data )
  k = 1; % we could still define how many functions should
         % be added to the space

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

  dummy    = [];

  enlarged = model.mexptr( 'rb_extension_PCA', model.mu, k );
end
%| \docupdate 
