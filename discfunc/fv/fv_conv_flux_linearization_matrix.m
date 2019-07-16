function flux_lin_mat = fv_conv_flux_linearization_matrix(model,model_data,U)
%function flux_lin_mat = fv_conv_flux_linearization_matrix(model,model_data,[U])
%
% function computing the linearized flux matrix of a convection problem. 
% simply reformatting the grid data suitably for pointwise evaluation
% by conv_flux_linearization. As evaluation points the points of suitable
% gauss-quadratures are chosen. The degree can be chosen in the
% model structure.
%
% required fields of model as required by conv_flux
%
% Function supports affine decomposition, i.e. different operation modes
% guided by optional field decomp_mode in model. See also the
% contents.txt for general explanation
%
% optional fields of model:
%   mu_names : names of fields to be regarded as parameters in vector mu
%   flux_quad_degree: degree of quadrature to be used for face-integrals.
%
% in 'coefficient' mode, U and the grid are empty.

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


% Bernard Haasdonk 11.4.2006

% determine affine_decomposition_mode as integer  
decomp_mode = model.decomp_mode;

grid = model_data.grid;

flux_lin_mat = [];

if decomp_mode > 0
  error('dont know how to handle decompositions in fv_conv_flux_linearization_matrix')
%  nfaces = size(grid.ECX,2);
%  UU = repmat(U,1,nfaces);
end;

if decomp_mode < 2
  nfaces = size(grid.ECX,2);
  UU = repmat(U,1,nfaces);
end;

if ~isempty(model.flux_quad_degree) 
  model.flux_quad_degree = 1;
end;

%if decomp_mode == 0
%    f = conv_flux(model, UU(:),grid.ECX(:), grid.ECY(:) );
% compute f by quadratures
[elids, edgeids] = ind2sub(size(grid.VI),1:length(grid.VI(:)));
PP = edge_quad_points(grid,elids, edgeids,model.flux_quad_degree);
[ff_lin, lambda] = ...
  model.conv_flux_derivative_ptr(PP, ...
           repmat(UU(:),model.flux_quad_degree,1), ...
           model);
f.Vx = edge_quad_eval_mean(grid,elids,edgeids,...
         model.flux_quad_degree,ff_lin(:,1));
f.Vy = edge_quad_eval_mean(grid,elids,edgeids,...
			   model.flux_quad_degree,ff_lin(:,2));
%f.lambda = ff_lin.lambda;
flux_lin_mat.Vx = reshape(f.Vx,size(UU));
flux_lin_mat.Vy = reshape(f.Vy,size(UU));
%elseif decomp_mode == 1
%  %    f = conv_flux(model, UU(:),grid.ECX(:), grid.ECY(:) );
%  [elids, edgeids] = ind2sub(size(grid.VI),1:length(grid.VI(:)));
%  PP = edge_quad_points(grid,elids, edgeids,model.flux_quad_degree);
%  ff_lin = conv_flux_linearization(model, repmat(UU(:),model.flux_quad_degree,1), ...
%				   PP(1,:), PP(2,:) );
%  f = cell(length(ff_lin),1);
%  for q = 1:length(ff);
%    f{q}.Fx = edge_quad_eval_mean(grid,elids,edgeids,...
%				  model.flux_quad_degree,ff{q}.Fx);
%    f{q}.Fy = edge_quad_eval_mean(grid,elids,edgeids,...
%				  model.flux_quad_degree,ff{q}.Fy);      
%    flux_lin_mat{q}.Fx = reshape(f{q}.Fx,size(UU));
%    flux_lin_mat{q}.Fy = reshape(f{q}.Fy,size(UU));
%  end;
%else % simply get coefficients from conv_flux and forward
%     %    flux_lin_mat = conv_flux(model, UU(:),grid.ECX(:), grid.ECY(:) );
%     flux_lin_mat = conv_flux(model, [],[],[] );
%end;

%if decomp_mode == 0
%  error('other modes than ''complete'' are currently not implemented! ');
%  flux_lin_mat.lambda = f.lambda;
%end;

 
