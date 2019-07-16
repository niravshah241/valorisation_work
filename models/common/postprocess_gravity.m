function data = postprocess_gravity(data, glob, params)
%function data = postprocess_gravity(data, glob, params)
% subtracts a previously added addent induced by gravitational effects.
%
% In case of a numerical scheme for the Richard's equation it is possible to
% model gravity by adding a defect growing with height to the initial data
% function. This way, the resulting concentration gradient diffuses a down-ward
% flow is created.
%
% required fields of params:
%   - 'gravity':   gravitational factor

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


data = data - glob(:,2) * params.gravity;

end
%| \docupdate 
