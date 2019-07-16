function OK = test_ldg_derivative
%function OK = test_ldg_derivative
%
% function testing, whether the ldg basis function derivatives
% gives approximately what is expected from finite differences.

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


% Bernard Haasdonk 28.8.2009

OK = 1;

for pdeg = 0:4
  params.nelements = 10;
  params.dimrange = 3;
  params.pdeg = pdeg;
%  df = ldgdiscfunc(params);
  
  for j = 1:100;
    x = rand(1,2);
    h = 1e-8;
    
    DPhis = ldg_evaluate_basis_derivative(x,params);
    
    % for each basis function compute finite differences
    Phis = ldg_evaluate_basis(x,params);
    Phis_xph = ldg_evaluate_basis(x+[h,0],params);
    Phis_yph = ldg_evaluate_basis(x+[0,h],params);
    
    nbasefunc = length(DPhis);
    for i = 1:nbasefunc;
      Dphi_appr = [Phis_xph(i,:)-Phis(i,:);  Phis_yph(i,:)-Phis(i,:)]'/h; 
      maxerr = max(max(abs(Dphi_appr-DPhis{i})));
      if maxerr>1e-3
	OK = 0;
	disp(['pdeg= ',num2str(pdeg),', i=',num2str(i),...
	      ',|DPhi-DPhi_appr|=',num2str(maxerr)])
	keyboard;
	error('ldg function derivatives not correct!');
      end;
    end;    
  end;
end;%| \docupdate 
