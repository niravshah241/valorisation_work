function varargout = rb_operators(model, detailed_data)
%function varargout = rb_operators(model, detailed_data)
%
% function calling the rb_operators method in the model
%
% Required fields of model:
%   rb_operators: handle
%

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


  [varargout{1:nargout}] = model.rb_operators(model, detailed_data);
end

