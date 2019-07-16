function res = power_vector2(x,pdeg)
%function res = power_vector2(x,pdeg)
%
% function computing the vector of all monomials of degree pdeg of
% the vector x, which is assumed to be 2-dimensional (resp: only
% first two entries are used)
% monomials up to deg 4 are explicitly implemented, hence fast,
% higher degrees are computed recursively. If using higher degree
% more often, simply insert explicit functions into case select list.

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


% Bernard Haasdonk 29.1.2009

switch pdeg
 case 0
  res =1;
 case 1
  res = [1;x(1);x(2)];
 case 2
  res = [1,x(1),x(2),x(1)^2, x(1)*x(2),x(2)^2]';
 case 3
  res = [1,x(1),x(2),x(1)^2, x(1)*x(2),x(2)^2,...
	 x(1)^3, x(1)^2*x(2), x(1) * x(2)^2, x(2)^3]';
 case 4
  res = [1,x(1),x(2),x(1)^2, x(1)*x(2),x(2)^2,...
	 x(1)^3, x(1)^2*x(2), x(1) * x(2)^2, x(2)^3,...
	 x(1)^4, x(1)^3*x(2), x(1)^2*x(2)^2, x(1)* x(2)^3, x(2)^4]';
 otherwise
  disp(['if using this pdeg more frequently, include in power_vector2', ...
       'explicitly for increased performance'])
  if pdeg~=uint8(pdeg)
    error('only integers as pdeg allowed!');
  end;  
  % recursive definition
  pv_pmin1 = power_vector2(x,pdeg-1);
  pv_new = x(1).^(pdeg:-1:0)'.* x(2).^(0:pdeg)';
  res = [pv_pmin1; pv_new];
end;
%| \docupdate 
