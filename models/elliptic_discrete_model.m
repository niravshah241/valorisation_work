function local_model = elliptic_discrete_model(model)
% function local_model = elliptic_discrete_model(model)
% 
% function creating a model with local functions out of a model
% with global functions, i.e. all volume-data functions can be
% evaluated in many elements simultanously by e.g.
%
% local_model.source(grid,eindices,loc,params)
%
% where eindices is a n times 1 vector of element-indices and loc
% is a local coordinate in the reference triangle. the result is a
% n times 1 vector of resulting values of the source.
% and boundary-data functions can be evaluated by e.g.
%
% local_model.dirichlet_values(grid,eindices,face_index,lloc,params)
%
% where eindices is a n times 1 vector of element-indices and lloc
% is a 1d coordinate in [0,1]. the result is a n times 1 vector of
% resulting dirichlet-values on the edge given by face_index.
%
% to avoid conflicts in the "global" model, all data functions can
% still be evaluated with (glob,params), where glob are global
% coordinates.

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


% Immanuel Maier, 25.03.2011

local_model = model;

% convert all data functions of model:

if model.has_reaction
  local_model.reaction = @(varargin) ...
      discrete_volume_values(model.reaction,varargin{:});
end;

if model.has_advection
  local_model.velocity = @(varargin) ...
      discrete_volume_values(model.velocity,varargin{:});
end;

if model.has_diffusivity
  local_model.diffusivity_tensor = @(varargin) ...
      discrete_volume_values(model.diffusivity_tensor,varargin{:});
end;

if model.has_source
 local_model.source = @(varargin) ...
     discrete_volume_values(model.source,varargin{:});
end;

if model.has_dirichlet_values
  local_model.dirichlet_values = @(varargin) ...
      discrete_boundary_values(model.dirichlet_values,varargin{:});
end;

if model.has_neumann_values
  local_model.neumann_values = @(varargin) ...
      discrete_boundary_values(model.neumann_values,varargin{:});
end;

if model.has_robin_values
  local_model.robin_alpha = @(varargin) ...
      discrete_boundary_values(model.robin_alpha,varargin{:});
  local_model.robin_beta = @(varargin) ...
      discrete_boundary_values(model.robin_beta,varargin{:});
  local_model.robin_values = @(varargin) ...
      discrete_boundary_values(model.robin_values,varargin{:});
end;

%%%%%%%%%%%% auxiliary functions:

function res = discrete_volume_values(func,varargin)

if nargin == 5
  loc = varargin{3};
  if numel(loc) == 0
    glob = [];
  else
    glob = local2global(varargin{1},varargin{2},loc);
  end;
elseif nargin == 3
  glob = varargin{1};
else
  error('incorrect number of inputs');
end;
res = func(glob,varargin{end});

function res = discrete_boundary_values(val,varargin)

if nargin == 6
  lloc = varargin{4};
  if numel(lloc) == 0
    glob = [];
  else
    loc = llocal2local(varargin{1},varargin{3},lloc);
    glob = local2global(varargin{1},varargin{2},loc);
  end;
elseif nargin == 3
  glob = varargin{1};
else
  error('incorrect number of inputs');
end;  
res = val(glob,varargin{end});
