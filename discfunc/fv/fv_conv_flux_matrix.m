function [flux_mat, lambda] = fv_conv_flux_matrix(model,model_data,U,conv_flux_ptr)
%function [flux_mat[, lambda]] = fv_conv_flux_matrix(model,model_data,[U,conv_flux_ptr])
%
% function computing the flux matrix of a convection problem. 
% simply reformatting the grid data suitably for pointwise evaluation
% by conv_flux. As evaluation points the points of suitable
% gauss-quadratures are chosen. The degree can be chosen in the
% model structure.
%
% in 'complete' (decomp_mode = 0) mode, flux_mat is a 
% 2 x nelements x nfaces matrix. In decomp_mode = 1 a cell array of
% such and in decomp_mode = 2 a simple vector of coefficients.
%
% Optionally, a divergence cleaning is performed. Ensure, that the
% divergence cleaning is not performed out of a neuman boundary, as these
% values will be reevaluated in setting of the boundary values, so an
% inconsistency would occur in that case. Instead lead the cleaning
% sweeps out of the domain by dirichlet boundaries.
% alternatively a divergence cleaning by optimization is possible.
%
% required fields of model as required by conv_flux
%
% optional fields of model:  
% divclean_mode: field, which indicates the divergence cleaning  
%                if it is set, further fields must be set:
%                possibility: 'none', 'sweep' or 'optimize'
%
% in case of 'sweep' additionally required:
% divclean_sweeps: cell array with information on divergence cleaning
%                  sweeps. Each element contains three fields:
%                  start_elements : vector of element indices used
%                                            as starting points for the sweep 
%                  steps          : number of rows/cols divergence cleansed 
%                  direction      : main direction of the block for cleaning
% 
% 
% Function supports affine decomposition, i.e. different operation modes
% guided by optional field decomp_mode in model. See also the
% contents.txt for general explanation

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

grid = [];
if ~isempty(model_data)
  grid = model_data.grid;
end
flux_mat = [];

if decomp_mode < 2
  nfaces = size(grid.ECX,2);
  UU = repmat(U,1,nfaces);
end;

if nargin < 4
  conv_flux_ptr = model.conv_flux_ptr;
end



%  if (~isequal(model.divclean_mode,'none')) & ...
%	isfield(model,'use_velocity_matrixfile') & ...
%	(model.use_velocity_matrixfile==1)
%    model.divclean_mode = 'none';
%    if model.verbose >= 5 
%      disp('divergence cleaning is skipped, as velocity from file is read!!');
%    end;
%  end;

if model.flux_linear && model.divclean_mode
  error('divergence clearning has to be adjusted!!')
%  % ensure that no affine decomposition is wanted!!!
%  if decomp_mode~=0
%    error('no simultaneous affine decomposition and divergence cleaning!');
%  end;
%  
%  % perform different cleaning methods:
%  f = conv_flux_linearization(model,UU(:),grid.ECX(:), grid.ECY(:) );  
%  V.Vx = reshape(f.Vx,size(grid.ECX));
%  V.Vy = reshape(f.Vy,size(grid.ECY));
%  
%  if isequal(model.divclean_mode,'sweep')
%    if model.verbose>9
%      disp('performing divergence cleaning by sweeping');
%    end;
%    Vclean = divergence_clean_sweep(model,model_data,V);
%  elseif isequal(model.divclean_mode,'optimize')
%    if model.verbose>9
%      disp('performing divergence cleaning by optimization');
%    end;
%    Vclean = divergence_clean_optimize(model, model_data,V); 
%  else
%    error('unknown divergence_cleaning mode!');  
%  end;
%  % set resulting flux values
%  flux_mat.Fx = Vclean.Vx.*UU;    
%  flux_mat.Fy = Vclean.Vy.*UU;    
%  %    flux_mat.Fx = reshape(Vclean.Vx, size(UU)).*UU;    
%  %    flux_mat.Fy = reshape(Vclean.Vy, size(UU)).*UU;    
else % no cleaning: ordinary flux values 
  if decomp_mode == 2 % simply get coefficients from conv_flux and forward
    flux_mat = conv_flux_ptr([],[],model);
  elseif decomp_mode == 0
    % compute f by quadratures
    [elids, edgeids] = ind2sub(size(grid.VI),1:length(grid.VI(:)));
    PP = edge_quad_points(grid,elids, edgeids, model.flux_quad_degree);
    %    ff = conv_flux(model,repmat(UU(:),model.flux_quad_degree,1), ...
    %             PP(1,:), PP(2,:));
    PP_nums = size(PP,1)/(size(grid.ECX,1)*size(grid.ECX,2));
    [ff, lambda] = conv_flux_ptr(PP,repmat(UU(:),PP_nums,1)', ...
                   model);
    fx = edge_quad_eval_mean(grid, elids, edgeids,...
                               model.flux_quad_degree, ff(:,1));
    fy = edge_quad_eval_mean(grid, elids, edgeids,...
                               model.flux_quad_degree, ff(:,2));
    f = [fx(:),fy(:)]';
% TODO: test if results in rows are faster
    flux_mat = reshape(f,[2,size(UU)]);
%    flux_mat.Fy = reshape(f.Fy,size(UU));
  elseif decomp_mode == 1
%    f = conv_flux(UU(:),grid.ECX(:), grid.ECY(:), model );
    [elids, edgeids] = ind2sub(size(grid.VI),1:length(grid.VI(:)));
    PP = edge_quad_points(grid,elids, edgeids,model.flux_quad_degree);
    PP_nums = size(PP,1)/(size(grid.ECX,1)*size(grid.ECX,2));
%    ff = conv_flux_ptr(model, repmat(UU(:),PP_nums,1), ...
%		   PP(1,:), PP(2,:));
    ff = conv_flux_ptr(PP, repmat(UU(:),PP_nums,1)',model);

    %f = cell(length(ff),1);
    for q = 1:length(ff);
      fx = edge_quad_eval_mean(grid,elids,edgeids,...
			       model.flux_quad_degree,ff{q}(:,1));
      fy = edge_quad_eval_mean(grid,elids,edgeids,...
			       model.flux_quad_degree,ff{q}(:,2));      
      f = [fx(:),fy(:)]';
      flux_mat{q} = reshape(f,[2,size(UU)]);
    end;
  end;
end;

% lambda only reasonable in case of no affine decomposition
%if decomp_mode == 0
%  flux_mat.lambda = f.lambda;
%end;
  

