function qres = demo_quadratures(pdeg,qdeg)
%function qres = demo_quadratures(pdeg,qdeg)
% small script demonstrating interval integration of polynomials
% of degree 'pdeg' with quadratures of degreed 'qdeg'
%
% Parameters:
%  pdeg: polynomial degree (default = 25)
%  qdeg: quadrature degree (default = 23)
%
% Return values:
%  qres: a matrix with the integration results

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

help demo_quadratures;

if nargin < 2
  qdeg = 23;
end;
if nargin < 1
  pdeg = 25;
end;
qres = zeros(pdeg+1,qdeg+1);
disp('results of integrating functions over unit interval:')
for d=0:pdeg
  f = @(x) x.^d;
  for q = 0:qdeg
    qres(d+1,q+1)=intervalquadrature(q,f);
    disp(['f(x) = x^',num2str(d),', qdeg = ',num2str(q),...
	   ' quadrature value = ',num2str(qres(d+1,q+1))]);
  end;
end;

