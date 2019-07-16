function U = fv_init_values(model,model_data)
%function U = fv_init_values(model,model_data)
%
% function computing the init-value DOF vector by evaluating the
% grid-center-of-gravity. In particular used for 1st order FV-schemes.
%
% required fields of model as used in model.init_values_ptr
%
% Function supports affine decomposition, i.e. different operation modes
% guided by optional field decomp_mode in params. See also the
% contents.txt for general explanation
%
% required fileds in model
%    init_values_ptr:   function pointer to init values function
%
% in 'coefficient' mode the model_data structure is empty

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


% Bernard Haasdonk 22.7.2006

if nargin ~= 2
  error('wrong number of parameters!');
end;

% affine decomposition is supported by init_data, so automatically
% satisfied.

decomp_mode = model.decomp_mode;

% TODO: The here called function init_values is still an old one that expects a
% "params" structure as last argument and uses string comparisons in order to
% choose the right initial data function. This part should therefore be
% rewritten with an implementation of function pointers.
if decomp_mode == 2
  U = model.init_values_ptr([],model);
else
  grid = model_data.grid;
  % initial values by midpoint evaluation
  U = model.init_values_ptr([grid.CX(:),grid.CY(:)],model);
end;

