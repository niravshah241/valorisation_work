function demo_edge_quad
%function demo_edge_quad
% demonstration of the edge-quadratures over a grid
% 
% a quartic divergence free velocity field is put over an
% unstructured triangular grid. Is is shown, how increasing the 
% quadrature degree improves the discrete divergence
%
% interestingly, a cubic velocity field is already integrated with
% 0 divergence by degree 2 Gaussian quadrature.

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


% Bernard Haadonk 13.5.2007

help demo_edge_quad

grid = triagrid;
%grid = rectgrid;

%params.xrange = [0,1];
%params.yrange = [0,1];
%params.xnumintervals = 5;
%params.ynumintervals = 5;
%grid = rectgrid(params);

params.axis_equal = 1;
params.no_lines = 1;

figure; 
%title('Comparison of edge quadrature errors');

% get all edges (even twice)
[elids, edgeids] = ind2sub(size(grid.VI),(1:length(grid.VI(:)))');
for degree = 1:3
  PP = edge_quad_points(grid,elids,edgeids,degree);
  VV = demo_velocity(PP);
  VX = edge_quad_eval(grid,elids,edgeids,degree,VV(:,1));
  VY = edge_quad_eval(grid,elids,edgeids,degree,VV(:,2));
  % compute discrete divergence
  VX = reshape(VX,size(grid.VI));
  VY = reshape(VY,size(grid.VI));
  D = sum(VX.*grid.NX + VY.*grid.NY,2);
%  Vav = sum(sqrt( (VX.*(grid.EL.^(-1))).^2 + ...
%		  (VY.*(grid.EL.^(-1))).^2), ...
%	    2)/grid.nneigh;
%  subplot(2,3,degree),plot_element_data(grid,Vav,params);
%  title(['deg = ',num2str(degree),...
%	 ', mean(|v|)']);  
  subplot(1,3,degree),plot_element_data(grid,D,params);
  title(['qdeg = ',num2str(degree),...
	 ', max |div v| = ', num2str(max(abs(D)))]);
%  keyboard;
end;

set(gcf,'Position',[28         297        1210         317]);

% demo function, evaluating a cubic velocity field in given points
% implements v(x,y)= (y^4,x^4) 
function FF = demo_velocity(PP)
FF = zeros(size(PP,1),2);
FF(:,1) = (2*PP(:,2)-1).^4;
FF(:,2) = (2*PP(:,1)-1).^4;

% TO BE ADJUSTED TO NEW SYNTAX
%| \docupdate 
