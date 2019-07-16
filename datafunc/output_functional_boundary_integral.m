function s = output_functional_boundary_integral(model,model_data,u)
%function s = output_functional_boundary_integral(model,model_data,u)
%
% function computing an output functional from the
% function U. Grid is needed to have the space discretization 
% information about U. A linear output function of the form of a
% boundary integral is computed:
%
%    output(u) = int_{boundary(Omega)} f u dx
%
% required fields of model:
%      output_function  : pointer to function f,
%                 e.g. constant 1 or 1/boundary_length, etc.
%                 computing the weight function f of the integral
%      output_integral_qdeg : degree of output integral computation
%
% u is assumed to be a function allowing a local evaluation
% u.evaluate(eindices, lcoord, grid,model)

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


% Bernard Haasdonk 27.2.2011

qdeg = model.output_integral_qdeg;
% find boundary edges and element indices
i = find(model_data.grid.NBI(:)<0)
[elind,edgeind] =  ind2sub(size(model_data.grid.NBI),i);
% define integral kernel, here x is a 1-d coordinate in [0,1] 
f = @(x) bnd_intkernel(x,elind,edgeind,model,model_data,u);
% compute output
res = intervalquadrature(qdeg,f);
s = sum(res);

%%%%%%%%%%%%%%% auxiliary function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% integral kernel for boundary integral
function res = bnd_intkernel(x,elind,edgeind,model,model_data,u)
res = [];
for local_edge_ind = 1:3;
  % determine edgelengths for transformation formula
  ind = find(edgeind==local_edge_ind);
  EL = model_data.grid.EL(elind(ind),local_edge_ind);
  lcoord = llocal2local(model_data.grid,local_edge_ind,x); % local coord for x on edge 1      
  glob = local2global(model_data.grid,elind(ind),lcoord);
  uval = u.evaluate(elind(ind), lcoord, model_data.grid, model);
  weightval = model.output_function(glob,model);
  res = [res;uval.*weightval.*EL]  
end;