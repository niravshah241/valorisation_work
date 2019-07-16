function s = output_functional_volume_integral(model,model_data,u)
%function s = output_functional_volume_integral(model,model_data,u)
%
% function computing an output functional from the discrete
% function U. Grid is needed to have the space discretization 
% information about U. A linear output function of the form of a
% volume integral is computed:
%
%    output(u) = int_Omega f u dx
%
% required fields of model:
%      output_function  : pointer to function f,
%                 e.g. output_function_box_mean.m
%                 computing the weight function f of the integral
%      output_integral_qdeg : degree of output integral computation
%
% u is assumed to be a function allowing a local evaluation

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


% Bernard Haasdonk 11.1.2011

qdeg = model.output_integral_qdeg;
eindices = 1:model_data.grid.nelements;
glob = local2global(model_data.grid, eindices, lcoord, model);

func = @(lcoord) ...
         u.evaluate(eindices, lcoord, model_data.grid, model).* ...
         model.output_function(glob);

local_int = triaquadrature(qdeg,func);

s = sum(local_int.*model_data.grid.A(eindices));


