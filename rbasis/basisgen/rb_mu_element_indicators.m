function [eta, eta_info] = rb_mu_element_indicators( detailed_data, ...
						  offline_data,...
						  MMesh, Delta_train, model)
%function [eta, eta_info] = rb_mu_element_indicators(
%                                detailed_data, [offline_data],
%                                MMesh, Delta_train, model)
%
% returns a vector with element-error indicators of the mu-grid for
% all leaf-elements. The current basis RB is used for estimator
% computation, the MMesh is the adaptive hierarchical mesh for the 
% parameters. Delta_train indicates the current maximum error
% estimator for the MMesh-vertices
%
% return value for a leaf element i is depending on the field 
% RB_element_indicator_mode in model. 
%   cog = center of gravity
%
% 'nodes':
%          eta(i) =  max_{mu in (V(e))}   Delta (Phi_t, mu)   
%
% 'nodes_cogs':
%          eta(i) =  max_{mu in (V(e) and cog(e))}  Delta (Phi_t, mu)   
%
% 'nodes_skippedrefs':
%          eta(i) =  max_{mu in (V(e))}   Delta (Phi_t, mu) +   
%                                     s(i)/s_max * max(Delta_train)
%
% 'nodes_cogs_skippedrefs':
%          eta(i) =  max_{mu in (V(e) and cog(e))}   Delta (Phi_t, mu) +   
%                                     s(i)/s_max * max(Delta_train)
%
% eta_info contains further detailed information, e.g. for debugging. 
% It consist of fields
%   Delta_train : vector of vertex indicators (i.e. Delta_train of input)
%   Delta_cog : vector of element-cog indicators in case of
%               nodes_cogs mode
%
% Required fields of model:
%     RB_element_indicator_mode : string indicating the mode, see above
%     RB_element_indicator_s_max : integer giving the maximum
%             acceptable number of skipped
%             refinement steps of elements.
%
% if Delta_train is empty, it is generated, if offline_data is
% empty, it is generated

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

% Bernard Haasdonk 22.3.2007

nleafelements = get(MMesh,'nleafelements');
eta = -1 * ones(1,nleafelements);
eta_info = [];

% determine all vertex error estimators

if isempty(offline_data)
%  model.Nmax = size(detailed_data.RB,2);
%  model.N = size(detailed_data.RB,2);
  offline_data = rb_offline_preparation(detailed_data,model);
end;

if isempty(Delta_train)
  M_train = get(MMesh,'vertex')';
  Delta_train = rb_test_indicator(detailed_data,...
				  offline_data,...
				  M_train,[],...
				  model);
end;

if isequal(model.RB_element_indicator_mode,'nodes_cogs') || ...
      isequal(model.RB_element_indicator_mode,'nodes_cogs_skippedrefs')
  
  % determine all cog vertex estimators
  M_cog = get(MMesh,'leafcogs')';
  Delta_cog = rb_test_indicator(model,detailed_data,...
				[],M_cog,[]);
end;

% assemble eta(i) =  max_{mu in (V(e) and cog(e)   Delta (Phi_t, mu)   

%keyboard;

li = get(MMesh,'leafelements');

lv = get(MMesh,'vertexindex',li);

Delta_vertex = reshape(Delta_train(lv(:)),size(lv));

Max_Delta_vertex = max(Delta_vertex,[],2);

% prepare return values
switch model.RB_element_indicator_mode
 case 'nodes_cogs'
  %          eta(i) =  max_{mu in (V(e) and cog(e))}  Delta (Phi_t, mu)    
  eta = max(Delta_cog(:),Max_Delta_vertex(:));
  eta_info.Delta_cog = Delta_cog;
 case 'nodes'
  %         eta(i) =  max_{mu in (V(e))}   Delta (Phi_t, mu)   
  %
  eta = Max_Delta_vertex(:);
 case 'nodes_skippedrefs'
  %          eta(i) =  max_{mu in (V(e))}   Delta (Phi_t, mu) +   
  %                                     s(i)/s_max * epsilon
  % determine number of skipped refinements of all leaf elements 
  s = MMesh.refine_steps - MMesh.creation_step(li);
  eta = Max_Delta_vertex(:) + ...
	max(Max_Delta_vertex) / ...
	model.RB_element_indicator_s_max * ...
	s(:);  
 case 'nodes_cogs_skippedrefs'
  %          eta(i) =  max_{mu in (V(e) and cogs)}   Delta (Phi_t, mu) +   
  %                                     s(i)/s_max * epsilon
  % determine number of skipped refinements of all leaf elements 
  s = MMesh.refine_steps - MMesh.creation_step(li);
  eta = max(Delta_cog(:),Max_Delta_vertex(:)) + ...
	max(Max_Delta_vertex) / ...
	model.RB_element_indicator_s_max * ...
	s(:);
  eta_info.Delta_cog = Delta_cog;
 otherwise 
  error('element_indicator_mode in grid-refinement unknown');
end;

eta_info.Delta_train = Delta_train;












