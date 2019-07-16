function [Gx, Gy] = evaluate_gradient(X,Y,p1data, params)
%function [Gx, Gy] = evaluate_gradient(X,Y,p1data, params)
%
% function performing an evaluation of the gradient of a given data
% field. Currently only cartesian grids are supported, i.e. the p1 function
% is bilinear. 
%
% required fields of params:
%
% xnumintervals: number of intervals in x-direction
% ynumintervals: number of intervals in y-direction
% xrange: interval in x-direction;
% yrange: interval in y-direction;

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

  
% Bernard Haasdonk 18.4.2006
 
  X = X(:);
  Y = Y(:);
  p1data = p1data(:);
  
  %Gx = zeros(size(X));  
  %Gy = zeros(size(X));  
  
  nx = params.xnumintervals;
  ny = params.ynumintervals;
  
  dx = (params.xrange(2)-params.xrange(1))/nx;
  dy = (params.yrange(2)-params.yrange(1))/ny;
  
  % determine the element indices in which the points reside (counting from
  % left lower domain row-wise)
  
  xi = floor((X-params.xrange(1))/dx)+1;
  yi = floor((Y-params.yrange(1))/dy)+1;
  el_ind = xi+(yi-1)*nx;
  % set left and upper border to last element
  i = find(abs(X-params.xrange(2))<eps);
  if ~isempty(i)
    xi(i) = nx;
  end;
  i = find(abs(Y-params.yrange(2))<eps);
  if ~isempty(i)
    yi(i) = ny;
  end;
  % check out of range
  if ~isempty(find( xi>nx | xi < 1 | yi > ny | yi < 1, 1 ))
    error('out of domain evaluation!');
  end;
  
  % determine the 4 vertex indices of the elements 
  % counting right lower corner 1 then counterclockwise
  v_ind = zeros(length(el_ind),4);
  v_ind(:,1) = (yi-1)*(nx+1)+1 +xi;
  v_ind(:,2) = (yi  )*(nx+1)+1 +xi;
  v_ind(:,3) = (yi  )*(nx+1)   +xi;
  v_ind(:,4) = (yi-1)*(nx+1)   +xi;
  
  % determine the fraction of coordinates ("modulo dx, dy")
  xfrac = X - params.xrange(1) - (xi-1)*dx;
  yfrac = Y - params.yrange(1) - (yi-1)*dy;
  
  % detAinv = (dx*dy)^-1
  % grad phi_1 = detAinv * [dy-y, -x   ]' 
  % grad phi_2 = detAinv * [   y,  x   ]' 
  % grad phi_3 = detAinv * [  -y,  dx-x]' 
  % grad phi_4 = detAinv * [y-dy,  x-dx]' 

  
  Gx = (dy-yfrac) .* p1data(v_ind(:,1)) ...
       + yfrac    .* p1data(v_ind(:,2)) ...
       - yfrac    .* p1data(v_ind(:,3)) ...
       +(yfrac-dy).* p1data(v_ind(:,4));
  
  Gy = - xfrac      .* p1data(v_ind(:,1)) ...
       + xfrac      .* p1data(v_ind(:,2)) ...
       +(dx - xfrac).* p1data(v_ind(:,3)) ...
       +(xfrac-dx)  .* p1data(v_ind(:,4));
  
  detAinv = (dx*dy)^-1;
  Gx = Gx * detAinv;
  Gy = Gy * detAinv;
  
  
% TO BE ADJUSTED TO NEW SYNTAX
%| \docupdate 
